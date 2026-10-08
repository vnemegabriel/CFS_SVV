#!/usr/bin/env python3
"""Copia un archivo reemplazando @nombre@ por el valor de la corrida, leído de run.ini.

    python3 render.py RUN_DIR ORIGEN DESTINO

Hace falta porque `-P` de Nektar++ cambia el parámetro nombrado pero no recalcula los
que dependen de él. Un @nombre@ sin valor es un error y la preparación falla.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg  # noqa: E402


def main(rd, src, dst):
    vals = dict(cfg.read_ini(os.path.join(rd, "run.ini"))["valores"])
    with open(src, encoding="utf8") as f:
        text = f.read()
    missing = sorted({m for m in re.findall(r"@(\w+)@", text) if m not in vals})
    if missing:
        sys.exit("render: %s sin valor para %s" % (src, ", ".join(missing)))
    with open(dst, "w", encoding="utf8") as f:
        f.write(re.sub(r"@(\w+)@", lambda m: vals[m.group(1)], text))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
