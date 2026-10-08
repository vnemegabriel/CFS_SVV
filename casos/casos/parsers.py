"""Salidas de cada solver: historial de fuerzas, tiempo de solver, celdas y firmas de falla.

Formatos verificados contra corridas reales (tests/fixtures):
  Nektar++ AeroForces  '# Time F1-press F1-visc F1-total F2-press ...'
  Nektar++ log         'Steps: N  Time: t  CPU Time: Xs'   tiempo de cada bloque de IO_InfoSteps
                       'Writing: "archivo" (Xs, XML)'       escritura, se descuenta
  OpenFOAM forces      '# Time total_x total_y total_z pressure_x ...'
  OpenFOAM log         'Time = t' ... 'ExecutionTime = X s  ClockTime = Y s'   acumulados
"""
import re

BAD_ROW = r"(?im)^\s*[-+]?\.?\d[^\n]*(?<![\w.])-?(?:nan|inf)(?!\w)"

PATTERNS = {
    "nektar": [
        ("divergio", r"NaN found"),
        ("divergio", BAD_ROW),
        ("fallo", r"Fatal\s*:|Segmentation fault|terminate called|MPI_ABORT|std::bad_alloc"),
    ],
    "openfoam": [
        ("divergio", r"sigFpe::sigHandler|(?i:residual = -?nan)"),
        ("divergio", BAD_ROW),
        ("fallo", r"FOAM FATAL|\[stack trace\]|Segmentation fault|MPI_ABORT|std::bad_alloc"),
    ],
}

COLUMNS = {"nektar": ("F1-total", "F2-total", "F3-total"),
           "openfoam": ("total_x", "total_y", "total_z")}


def forces(solver, paths):
    """Concatena archivos en orden; descarta filas incompletas y tiempos repetidos."""
    names = COLUMNS[solver]
    out = {"t": [], "fx": [], "fy": [], "fz": []}
    for path in paths:
        header = None
        with open(path, encoding="utf8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if s.startswith("#"):
                    tokens = s.lstrip("#").split()
                    if tokens[:1] == ["Time"]:
                        header = tokens
                        if names[0] not in header:
                            raise ValueError("%s: no tiene la columna %s" % (path, names[0]))
                    continue
                parts = s.replace("(", " ").replace(")", " ").split()
                if header is None or len(parts) != len(header):
                    continue
                try:
                    row = dict(zip(header, map(float, parts)))
                except ValueError:
                    continue
                if out["t"] and not row["Time"] > out["t"][-1]:
                    continue
                out["t"].append(row["Time"])
                for key, name in zip(("fx", "fy", "fz"), names):
                    out[key].append(row.get(name, 0.0))
    return out


def timing(solver, path):
    """(tiempos, segundos de solver acumulados sin escritura a disco)."""
    with open(path, encoding="utf8", errors="replace") as f:
        lines = f.readlines()
    return _nektar(lines) if solver == "nektar" else _openfoam(lines)


def _nektar(lines):
    step = re.compile(r"Steps:\s*\d+\s+Time:\s*(\S+)\s+CPU Time:\s*(\S+?)s\b")
    write = re.compile(r'Writing: ".*?" \((\S+?)s')
    t, sec, total, io = [0.0], [0.0], 0.0, 0.0
    for line in lines:
        m = write.search(line)
        if m:
            io += float(m.group(1))
            continue
        m = step.search(line)
        if m:
            total += max(0.0, float(m.group(2)) - io)
            io = 0.0
            t.append(float(m.group(1)))
            sec.append(total)
    return (t, sec) if len(t) > 1 else ([], [])


def _openfoam(lines):
    clock = re.compile(r"^ExecutionTime = (\S+) s\s+ClockTime = (\S+) s")
    now = re.compile(r"^Time = (\S+)")
    t, sec, start, current, looping = [], [], 0.0, None, False
    for line in lines:
        if line.startswith("Starting time loop"):
            looping = True
            continue
        m = now.match(line)
        if m:
            current = float(m.group(1).rstrip("s"))
            continue
        m = clock.match(line)
        if m:
            wall = float(m.group(2))
            if not looping or current is None:
                start = wall
            else:
                t.append(current)
                sec.append(wall - start)
    return t, sec


def errors(path):
    """Normas contra ExactSolution que Nektar++ imprime al terminar: {'L2_rho': v, 'Linf_rho': v}."""
    rx = re.compile(r"^L (2|inf) error \(variable (\w+)\)\s*:\s*(\S+)")
    out = {}
    with open(path, encoding="utf8", errors="replace") as f:
        for line in f:
            m = rx.match(line.strip())
            if m:
                out["L%s_%s" % (m.group(1), m.group(2))] = float(m.group(3))
    return out


def count_cells(solver, path):
    try:
        if solver == "nektar":
            with open(path, encoding="utf8", errors="replace") as f:
                block = re.search(r"<ELEMENT[^>]*>(.*?)</ELEMENT>", f.read(), re.S)
            if not block or "COMPRESSED" in block.group(0)[:200]:
                return None
            return len(re.findall(r"<[QTHRAP] ID=", block.group(1))) or None
        owner = path.rstrip("/") + "/constant/polyMesh/owner"
        with open(owner, encoding="utf8", errors="replace") as f:
            m = re.search(r"nCells:\s*(\d+)", f.read(4000))
        return int(m.group(1)) if m else None
    except OSError:
        return None
