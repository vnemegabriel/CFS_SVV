#!/usr/bin/env python3
"""Casos de simulación: planificar, verificar, correr, vigilar y post-procesar.

  casos.py plan    CASO [--force]           matrix.txt desde [sweep]
  casos.py check   CASO                     verifica todo sin correr nada
  casos.py run     CASO [--solo PATRÓN] [--retry]
  casos.py submit  CASO [--solo PATRÓN] [--retry]
  casos.py exec    CASO ID                  una corrida; la usan los trabajos en cola
  casos.py status  CASO [--todas]
  casos.py post    CASO [--desde RAÍZ ...] [--reprocesar] [--cada SEGUNDOS]
  casos.py unlock  CASO ID

--maquina NOMBRE elige la sección de maquinas.ini; sin eso se busca por host.
"""
import argparse
import os
import shlex
import signal
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cfg       # noqa: E402
import parsers   # noqa: E402
import post      # noqa: E402
import report    # noqa: E402
import runner    # noqa: E402
import telemetry  # noqa: E402


def interrupt(signum, frame):
    raise KeyboardInterrupt


def refresh(caso, roots):
    try:
        report.regenerate(caso, roots)
    except Exception as e:
        print("!! reporte: %s: %s" % (type(e).__name__, e), flush=True)


def pending(caso, mach, args):
    for row in caso.rows(args.solo):
        state = runner.read_status(caso.run_dir(row["id"], mach))["estado"]
        if state in ("ok", "corriendo", "encolada"):
            continue
        if state in runner.TERMINAL and not args.retry:
            continue
        yield row


# ---------------------------------------------------------------- comandos

def cmd_plan(caso, args):
    rows = cfg.plan(caso)
    if not rows:
        sys.exit("error: [sweep] vacío")
    if os.path.exists(caso.matrix_path) and not args.force:
        sys.exit("error: %s ya existe (--force para rehacerla)" % caso.matrix_path)
    cfg.write_matrix(caso.matrix_path, rows)
    print("%d corridas -> %s" % (len(rows), caso.matrix_path))


def cmd_check(caso, mach, args):
    errors, warnings = [], []
    if caso.solver not in parsers.PATTERNS:
        errors.append("[caso] solver '%s': se espera %s" % (caso.solver, " o ".join(parsers.PATTERNS)))
    for sec in ("caso", "values", "run", "post"):
        if not caso.ini.has_section(sec):
            errors.append("falta la sección [%s]" % sec)
    try:
        rows = caso.rows()
    except (OSError, ValueError) as e:
        errors.append(str(e))
        rows = []

    seen, first = {}, None
    for row in rows:
        rd = caso.run_dir(row["id"], mach)
        try:
            vals = cfg.resolve(caso, row, mach)
            first = first or (vals, rd)
            for sec, key in (("run", "prepare"), ("run", "command"), ("run", "watch"),
                             ("post", "forces"), ("post", "log"), ("post", "mesh"),
                             ("post", "unit"), ("post", "window"), ("post", "tol_cd"), ("post", "tol_cl")):
                cfg.template(caso, key, vals, rd, sec)
            if post.kind(caso) == "forces":
                missing = [k for k in post.REQUIRED if not post.isnum(vals.get(k))]
                if missing:
                    raise ValueError("[values] sin valor numérico para " + ", ".join(missing))
                try:
                    float(cfg.template(caso, "window", vals, rd, "post") or "")
                except ValueError:
                    raise ValueError("[post] window debe ser un número")
                if not any(cfg.template(caso, "tol_" + q, vals, rd, "post") for q in ("cd", "cl")):
                    raise ValueError("[post] necesita tol_cd o tol_cl")
            mesh = cfg.template(caso, "mesh", vals, rd, "post")
            if mesh and not mesh.startswith(rd) and not os.path.exists(mesh):
                warnings.append("%s: no existe la malla %s" % (row["id"], mesh))
        except ValueError as e:
            seen.setdefault(str(e), []).append(row["id"])
    for msg, ids in seen.items():
        errors.append("%s  (%d corridas, p.ej. %s)" % (msg, len(ids), ids[0]))

    if first and caso.get("run", "probe"):
        probe = cfg.template(caso, "probe", *first)
        r = subprocess.run(runner.shell(probe, mach, caso.solver), capture_output=True, text=True)
        if r.returncode != 0:
            errors.append("[run] probe falló en esta máquina: %s" % (r.stderr.strip() or probe))
    if mach["taubench_s"] is None:
        warnings.append("máquina %s sin taubench_s: el costo normalizado queda vacío hasta calibrar"
                        % mach["nombre"])
    if caso.get("report", "kind") not in report.KINDS:
        warnings.append("[report] kind debería ser " + " o ".join(report.KINDS))
    if post.kind(caso) not in ("forces", "errores"):
        errors.append("[post] kind debe ser forces o errores")
    if caso.ini.has_section("telemetry"):
        flight, note = telemetry.flight_cd(dict(caso.items("telemetry")), caso.dir)
        (warnings if not flight else []).append(note)

    print("caso %s   máquina %s   %d corridas" % (caso.name, mach["nombre"], len(rows)))
    for e in warnings:
        print("  aviso  " + e)
    for e in errors:
        print("  ERROR  " + e)
    print("ok" if not errors else "%d errores" % len(errors))
    return 1 if errors else 0


def cmd_run(caso, mach, args):
    todo = list(pending(caso, mach, args))
    print("%d corridas por correr" % len(todo), flush=True)
    for row in todo:
        print("== %s" % row["id"], flush=True)
        state = runner.execute(caso, row, mach)
        print("   %s  %s" % (state, runner.read_status(caso.run_dir(row["id"], mach))["motivo"]), flush=True)
        refresh(caso, [mach["runs_root"]])


def cmd_submit(caso, mach, args):
    submit = mach.get("submit", "local")
    if submit == "local":
        sys.exit("error: la máquina %s no tiene submit; usar `run`" % mach["nombre"])
    header = "#!/bin/bash"
    if mach.get("job_header"):
        with open(os.path.expanduser(mach["job_header"]), encoding="utf8") as f:
            header = f.read()
    for row in pending(caso, mach, args):
        rd = caso.run_dir(row["id"], mach)
        os.makedirs(rd, exist_ok=True)
        procs = cfg.number(row.get("procs", mach["procs"]))
        fill = lambda text: (text.replace("{id}", row["id"]).replace("{procs}", str(procs))
                             .replace("{run}", rd))
        job = os.path.join(rd, "job.sh")
        with open(job, "w", encoding="utf8") as f:
            f.write(fill(header).rstrip() + "\n")
            f.write("exec python3 %s --maquina %s exec %s %s\n" % (
                shlex.quote(os.path.join(cfg.HERE, "casos.py")), shlex.quote(mach["nombre"]),
                shlex.quote(caso.dir), shlex.quote(row["id"])))
        r = subprocess.run(["bash", "-c", submit.replace("{job}", shlex.quote(job))],
                           cwd=rd, capture_output=True, text=True)
        if r.returncode == 0:
            runner.write_status(rd, "encolada", r.stdout.strip())
            print("encolada  %s  %s" % (row["id"], r.stdout.strip()))
        else:
            print("!! %s: envío falló: %s" % (row["id"], r.stderr.strip()))


def cmd_exec(caso, mach, args):
    rows = [r for r in caso.rows() if r["id"] == args.id]
    if not rows:
        sys.exit("error: %s no está en matrix.txt" % args.id)
    print(runner.execute(caso, rows[0], mach), flush=True)
    refresh(caso, [mach["runs_root"]])


def cmd_status(caso, mach, args):
    rows = report.gather(caso, [mach["runs_root"]])
    count = {}
    for r in rows:
        count[r["estado"]] = count.get(r["estado"], 0) + 1
        if args.todas or r["estado"] != "pendiente":
            print("%-34s %-14s %s" % (r["id"], r["estado"], r["motivo"][:100]))
    print("   ".join("%s %d" % kv for kv in sorted(count.items())))


def cmd_post(caso, mach, args):
    roots = [os.path.expanduser(r) for r in (args.desde or [mach["runs_root"]])]
    while True:
        try:
            for row in caso.rows():
                for root in roots:
                    rd = os.path.join(root, caso.name, row["id"])
                    if not os.path.exists(os.path.join(rd, "run.ini")):
                        continue
                    if runner.read_status(rd)["estado"] not in runner.TERMINAL:
                        continue
                    if args.reprocesar or not os.path.exists(os.path.join(rd, "result.ini")):
                        runner.repost(caso, rd)
            print(report.regenerate(caso, roots), flush=True)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print("!! %s: %s" % (type(e).__name__, e), flush=True)
        if not args.cada:
            return
        time.sleep(args.cada)


def cmd_unlock(caso, mach, args):
    rd = caso.run_dir(args.id, mach)
    runner.release(rd)
    if runner.read_status(rd)["estado"] in ("corriendo", "encolada"):
        runner.write_status(rd, "interrumpida", "desbloqueada a mano")
    print("desbloqueada %s" % args.id)


# ---------------------------------------------------------------- entrada

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maquina")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "check", "run", "submit", "exec", "status", "post", "unlock"):
        p = sub.add_parser(name)
        p.add_argument("caso", metavar="CASO")
        if name in ("exec", "unlock"):
            p.add_argument("id")
        if name in ("run", "submit"):
            p.add_argument("--solo", metavar="PATRÓN")
            p.add_argument("--retry", action="store_true", help="repetir las que no terminaron ok")
        if name == "plan":
            p.add_argument("--force", action="store_true")
        if name == "status":
            p.add_argument("--todas", action="store_true")
        if name == "post":
            p.add_argument("--desde", nargs="+", metavar="RAÍZ")
            p.add_argument("--reprocesar", action="store_true")
            p.add_argument("--cada", type=float, metavar="SEGUNDOS")
    args = ap.parse_args(argv)

    signal.signal(signal.SIGTERM, interrupt)
    signal.signal(signal.SIGHUP, interrupt)
    try:
        caso = cfg.Caso(args.caso)
        if args.cmd == "plan":
            return cmd_plan(caso, args)
        mach = cfg.machine(args.maquina)
        return globals()["cmd_" + args.cmd](caso, mach, args) or 0
    except KeyboardInterrupt:
        print("\ninterrumpido: la corrida en curso quedó como 'interrumpida'", file=sys.stderr)
        return 130
    except (LookupError, OSError, ValueError) as e:
        print("error: %s" % e, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
