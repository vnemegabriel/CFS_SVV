#!/usr/bin/env python3
"""View a Nektar++ VTU time series (.pvd) produced by postNek.py.

  python viewNek.py vtu/naca.pvd --field u              # interactive
  python viewNek.py vtu/naca.pvd --field u --movie w.mp4
  python viewNek.py vtu/naca.pvd --field u --shot t.png --zoom 4
"""
import argparse
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import pyvista as pv


def series(pvd: Path) -> list:
    root = ET.parse(pvd).getroot()
    sets = root.findall(".//DataSet")
    return [(float(d.get("timestep")), pvd.parent / d.get("file")) for d in sets]


def make_plotter(off_screen: bool) -> pv.Plotter:
    pl = pv.Plotter(off_screen=off_screen)
    pl.set_background("white")
    return pl


def draw(pl: pv.Plotter, mesh, field: str, clim) -> None:
    pl.add_mesh(mesh, scalars=field, cmap="turbo", clim=clim,
                smooth_shading=True, show_edges=False,
                scalar_bar_args={"color": "black", "title": field})


def limits(frames: list, field: str) -> tuple:
    lo = min(f.get_array(field).min() for f in frames)
    hi = max(f.get_array(field).max() for f in frames)
    return float(lo), float(hi)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("pvd")
    p.add_argument("--field", default="u")
    p.add_argument("--movie", default=None, help="write mp4/gif instead of showing")
    p.add_argument("--shot", default=None, help="write a png of the last frame")
    p.add_argument("--zoom", type=float, default=1.0)
    a = p.parse_args()

    steps = series(Path(a.pvd).resolve())
    frames = [pv.read(f) for _, f in steps]
    if a.field not in frames[0].array_names:
        print(f"fields: {frames[0].array_names}", file=sys.stderr)
        return 1
    clim = limits(frames, a.field)

    pl = make_plotter(off_screen=bool(a.movie or a.shot))
    draw(pl, frames[-1], a.field, clim)
    pl.view_xy()
    pl.camera.zoom(a.zoom)

    if a.shot:
        pl.screenshot(a.shot)
        return 0
    if not a.movie:
        pl.show()
        return 0

    open_movie = pl.open_gif if a.movie.endswith(".gif") else pl.open_movie
    open_movie(a.movie)
    for (t, _), mesh in zip(steps, frames):
        pl.clear()
        draw(pl, mesh, a.field, clim)
        pl.add_text(f"t = {t:.3f}", color="black", font_size=10)
        pl.write_frame()
    pl.close()
    print(f"wrote {a.movie}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
