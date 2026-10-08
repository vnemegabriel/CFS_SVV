#!/usr/bin/env python3
"""Batch FieldConvert for Nektar++ cases -> VTU time series (.pvd).

Runs on Windows (through WSL) and on Linux (native) with the same command.

  python3 nek_post.py --mesh naca.xml --session session_naca_ex-Complete.xml
  python3 nek_post.py ... --order 4          # high-order VTU (Lagrange cells)
  python3 nek_post.py ... --modules mach pressure vorticity
"""
import argparse
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

WSL = sys.platform.startswith("win")


def to_exec_path(p: Path) -> str:
    """Path as seen by the process that runs FieldConvert."""
    if not WSL:
        return str(p)
    out = subprocess.run(["wsl", "wslpath", "-a", str(p).replace("\\", "/")],
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def run(args: list) -> None:
    cmd = ["wsl", "-e"] + args if WSL else args
    if subprocess.run(cmd).returncode != 0:
        raise SystemExit("FieldConvert failed:\n  " + " ".join(cmd))


def require(*paths: Path) -> None:
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise SystemExit("not found:\n  " + "\n  ".join(missing))


def chk_time(chk: Path, fallback: float) -> float:
    info = chk / "Info.xml" if chk.is_dir() else None
    if info is None or not info.exists():
        return fallback
    node = ET.parse(info).getroot().find(".//Time")
    return float(node.text) if node is not None else fallback


def chk_index(chk: Path) -> int:
    m = re.search(r"_(\d+)\.chk$", chk.name)
    return int(m.group(1)) if m else 10**9


def collect(case: Path, stem: str) -> list:
    """[(time, path)] over checkpoints plus the final field, duplicates dropped."""
    srcs = sorted(case.glob(f"{stem}_*.chk"), key=chk_index)
    final = case / f"{stem}.fld"
    srcs += [final] if final.exists() else []
    steps, seen = [], set()
    for i, src in enumerate(srcs):
        t = chk_time(src, float(i))
        if t not in seen:
            seen.add(t)
            steps.append((t, src))
    return steps


def convert(case: Path, mesh: Path, session: Path, src: Path, dst: Path,
            modules: list, order: int, npts: int) -> None:
    spec = to_exec_path(dst)
    if order:
        spec += ":vtu:highorder"
    cmd = ["FieldConvert", "-f"]
    for m in modules:
        cmd += ["-m", m]
    if npts:
        cmd += ["-n", str(npts)]
    cmd += [to_exec_path(mesh), to_exec_path(session), to_exec_path(src), spec]
    run(cmd)


def write_pvd(pvd: Path, entries: list) -> None:
    rows = "\n".join(
        f'    <DataSet timestep="{t:.10g}" group="" part="0" file="{f}"/>'
        for t, f in entries)
    pvd.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        f'  <Collection>\n{rows}\n  </Collection>\n</VTKFile>\n')


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--case", default=".", help="case directory")
    p.add_argument("--mesh", required=True, help="geometry xml, e.g. naca.xml")
    p.add_argument("--session", required=True, help="session xml")
    p.add_argument("--stem", default=None, help="field basename (default: mesh stem)")
    p.add_argument("--out", default="vtu", help="output subdirectory")
    p.add_argument("--modules", nargs="*", default=[], help="FieldConvert -m modules")
    p.add_argument("--order", type=int, default=0, help="write high-order VTU")
    p.add_argument("--npts", type=int, default=0, help="equispaced points per direction")
    a = p.parse_args()

    case = Path(a.case).resolve()
    mesh, session = case / a.mesh, case / a.session
    stem = a.stem or mesh.stem
    outdir = case / a.out
    outdir.mkdir(exist_ok=True)

    require(mesh, session)
    sources = collect(case, stem)
    if not sources:
        print(f"no {stem}_*.chk or {stem}.fld in {case}", file=sys.stderr)
        return 1

    entries = []
    for i, (t, src) in enumerate(sources):
        dst = outdir / f"{stem}_{i:04d}.vtu"
        convert(case, mesh, session, src, dst, a.modules, a.order, a.npts)
        entries.append((t, dst.name))

    write_pvd(outdir / f"{stem}.pvd", entries)
    print(f"{len(entries)} steps -> {outdir / (stem + '.pvd')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
