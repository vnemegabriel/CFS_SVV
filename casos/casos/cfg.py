"""Configuración: máquina, caso, matriz y valores de una corrida.

Todo es texto plano que una persona puede leer y editar:
  maquinas.ini              una sección por máquina
  <caso>/caso.ini    qué correr y cómo post-procesarlo
  <caso>/matrix.txt      una corrida por línea, generada por `casos.py plan`
"""
import configparser
import fnmatch
import itertools
import math
import os
import socket

HERE = os.path.dirname(os.path.abspath(__file__))

SAFE = {k: getattr(math, k) for k in
        ("sqrt", "exp", "log", "log10", "sin", "cos", "tan", "atan2", "radians", "degrees", "pi")}
SAFE.update(abs=abs, min=min, max=max, round=round, int=int, float=float, str=str)


def read_ini(path):
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=("#",))
    cp.optionxform = str
    with open(path, encoding="utf8") as f:
        cp.read_file(f)
    return cp


def read_ini_flat(path):
    if not os.path.exists(path):
        return {}
    ini = read_ini(path)
    return {"%s.%s" % (s, k): v for s in ini.sections() for k, v in ini[s].items()}


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return "%.10g" % v
    return " ".join(str(v).split()) if "\n" in str(v) else str(v)


def write_ini(path, sections):
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    for sec, kv in sections.items():
        cp[sec] = {k: fmt(v) for k, v in kv.items()}
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf8") as f:
        cp.write(f)
    os.replace(tmp, path)


def number(s):
    for cast in (int, float):
        try:
            return cast(s)
        except (TypeError, ValueError):
            pass
    return s


def evaluate(expr, env):
    return eval(expr, {"__builtins__": {}}, env)


# ---------------------------------------------------------------- máquina

def machine(name=None):
    ini = read_ini(os.environ.get("CAMP_MAQUINAS") or os.path.join(HERE, "maquinas.ini"))
    name = name or os.environ.get("CAMP_MAQUINA")
    host = socket.gethostname()
    for sec in ini.sections():
        hosts = ini[sec].get("hosts", sec).split()
        if name == sec or (not name and any(fnmatch.fnmatch(host, h) for h in hosts)):
            m = dict(ini[sec])
            m["nombre"] = sec
            m["procs"] = int(m["procs"])
            m["runs_root"] = os.path.expanduser(m["runs_root"])
            m["taubench_s"] = float(m["taubench_s"]) if m.get("taubench_s") else None
            m["vars"] = {k[4:]: v for k, v in m.items() if k.startswith("var_")}
            return m
    raise LookupError("maquinas.ini no tiene una sección para %s" % (name or "el host " + host))


def taubench_of(name):
    try:
        return machine(name)["taubench_s"]
    except (LookupError, OSError):
        return None


# ---------------------------------------------------------------- caso

class Caso:
    def __init__(self, name):
        self.dir = os.path.abspath(name if os.path.isdir(name) else os.path.join(HERE, name))
        path = os.path.join(self.dir, "caso.ini")
        if not os.path.exists(path):
            raise FileNotFoundError("no existe %s" % path)
        self.ini = read_ini(path)
        self.name = os.path.basename(self.dir)
        self.solver = self.get("caso", "solver")
        self.matrix_path = os.path.join(self.dir, "matrix.txt")
        self.results_dir = os.path.join(self.dir, "resultados")

    def get(self, sec, key, default=None):
        if self.ini.has_option(sec, key) and self.ini[sec][key].strip():
            return self.ini[sec][key].strip()
        return default

    def items(self, sec):
        return list(self.ini[sec].items()) if self.ini.has_section(sec) else []

    def run_dir(self, run_id, mach, root=None):
        return os.path.join(root or mach["runs_root"], self.name, run_id)

    def rows(self, pattern=None):
        rows = read_matrix(self.matrix_path)
        return [r for r in rows if not pattern or fnmatch.fnmatch(r["id"], pattern)]


def plan(caso):
    """Producto de las claves de [sweep]. `mach,alt = 1.8,3000 1.2,5200` son columnas que van juntas."""
    groups = []
    for key, value in caso.items("sweep"):
        names = key.split(",")
        combos = [tuple(v.split(",")) for v in value.split()]
        for c in combos:
            if len(c) != len(names):
                raise ValueError("[sweep] %s: '%s' no tiene %d valores" % (key, ",".join(c), len(names)))
        groups.append([dict(zip(names, c)) for c in combos])
    pattern = caso.get("caso", "id")
    rows = []
    for combo in itertools.product(*groups):
        row = {}
        for part in combo:
            row.update(part)
        rows.append(dict(id=pattern.format(**row), **row))
    return rows


def write_matrix(path, rows):
    cols = list(rows[0])
    width = {c: max([len(c)] + [len(r[c]) for r in rows]) for c in cols}
    line = lambda r: "  ".join(r[c].ljust(width[c]) for c in cols).rstrip() + "\n"
    with open(path, "w", encoding="utf8") as f:
        f.write("# Una corrida por línea. Borrar una línea la saca del caso.\n")
        f.write(line({c: c for c in cols}))
        for r in rows:
            f.write(line(r))


def read_matrix(path):
    if not os.path.exists(path):
        raise FileNotFoundError("no existe %s: correr `casos.py plan` primero" % path)
    header, rows = None, []
    with open(path, encoding="utf8") as f:
        for n, raw in enumerate(f, 1):
            parts = raw.split("#", 1)[0].split()
            if not parts:
                continue
            if header is None:
                header = parts
                continue
            if len(parts) != len(header):
                raise ValueError("matrix.txt línea %d: %d columnas, se esperaban %d"
                                 % (n, len(parts), len(header)))
            rows.append(dict(zip(header, parts)))
    ids = [r["id"] for r in rows]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise ValueError("matrix.txt: ids repetidos: " + ", ".join(dup))
    return rows


# ---------------------------------------------------------------- valores de una corrida

def resolve(caso, row, mach):
    """Máquina, después fila de la matriz, después [values] en orden. La fila manda."""
    pending = [k for k, v in row.items() if v == "PENDIENTE"]
    if pending:
        raise ValueError("matrix.txt: %s = PENDIENTE" % ", ".join(pending))
    vals = {"procs": mach["procs"]}
    vals.update(mach.get("vars", {}))
    vals.update({k: number(v) for k, v in row.items()})
    vals["id"] = row["id"]
    env = dict(SAFE, **vals)
    for key, expr in caso.items("values"):
        if key in row:
            continue
        try:
            vals[key] = env[key] = evaluate(expr, env)
        except Exception as e:
            raise ValueError("[values] %s = %s  ->  %s: %s" % (key, expr, type(e).__name__, e))
    return vals


def template(caso, key, vals, rd, section="run"):
    text = caso.get(section, key)
    if not text:
        return None
    ctx = dict(vals, run=rd, case=os.path.join(rd, "case"), carpeta=caso.dir, tools=HERE)
    try:
        return text.format_map(ctx)
    except KeyError as e:
        raise ValueError("[%s] %s: no hay valor para {%s}" % (section, key, e.args[0]))
