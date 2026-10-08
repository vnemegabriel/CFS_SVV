"""Salida de una corrida -> result.ini normalizado.

Dos tipos, según [post] kind:
  forces   historial de fuerzas -> punto de parada -> Cd, Cl (validación, aplicación, costo)
  errores  normas contra ExactSolution al final de la corrida (verificación, calibración)

Normalizado significa: mismos coeficientes calculados con la misma referencia
(rho, U, aref, lref de [values]) sin importar el solver, tiempo en unidades
convectivas y costo dividido por el tiempo de TauBench de la máquina.
"""
import glob
import math
import os

import cfg
import num
import parsers

REQUIRED = ("rho", "U", "aref", "lref")


def isnum(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def kind(caso):
    return caso.get("post", "kind", "forces")


def coefficients(forces, vals):
    a, b = math.radians(vals.get("alpha", 0.0)), math.radians(vals.get("beta", 0.0))
    drag = (math.cos(a) * math.cos(b), math.sin(a) * math.cos(b), math.sin(b))
    lift = (-math.sin(a), math.cos(a), 0.0)
    qa = 0.5 * vals["rho"] * vals["U"] ** 2 * vals["aref"]
    f = list(zip(forces["fx"], forces["fy"], forces["fz"]))
    return {"cd": [sum(c * d for c, d in zip(v, drag)) / qa for v in f],
            "cl": [sum(c * d for c, d in zip(v, lift)) / qa for v in f]}


def process(caso, rd, vals, meta):
    path = lambda key: cfg.template(caso, key, vals, rd, "post")
    solver = caso.solver
    res = {"corrida": {k: meta.get(k, "") for k in ("id", "caso", "maquina", "host", "fecha")}}
    res["corrida"].update(solver=solver, procesos=vals["procs"])
    res["malla"] = mesh(caso, vals, path("mesh"))
    logs = sorted(glob.glob(path("log") or ""))

    if kind(caso) == "errores":
        res["errores"] = parsers.errors(logs[0]) if logs else {}
        res["costo"] = cost(solver, logs, None, vals["procs"], meta)
        return res

    missing = [k for k in REQUIRED if not isnum(vals.get(k))]
    if missing:
        raise ValueError("faltan en [values] para normalizar: " + ", ".join(missing))
    res.update(history(caso, vals, sorted(glob.glob(path("forces") or "")), path))
    res["costo"] = cost(solver, logs, res["parada"].get("t_parada"), vals["procs"], meta)
    return res


def mesh(caso, vals, meshpath):
    out = {}
    if not meshpath:
        return out
    cells = parsers.count_cells(caso.solver, meshpath)
    out["celdas"] = cells
    expr = caso.get("post", "dof")
    if cells and expr:
        dof = cfg.evaluate(expr, dict(cfg.SAFE, celdas=cells, **vals))
        out.update(dof=dof, h_eff=dof ** (-1.0 / int(caso.get("post", "dim", "2"))))
    return out


def history(caso, vals, files, path):
    forces = parsers.forces(caso.solver, files) if files else None
    n = len(forces["t"]) if forces else 0
    out = {"historia": {"filas": n, "archivo": " ".join(files)}, "parada": {}, "resultado": {}}
    if not n:
        out["parada"] = {"cumplio": "no", "nota": "sin historial de fuerzas"}
        return out

    series = coefficients(forces, vals)
    t = forces["t"]
    unit = path("unit") or "conv"
    x = [ti * vals["U"] / vals["lref"] for ti in t] if unit == "conv" else list(t)
    window = float(path("window"))
    tols = {q: float(path("tol_" + q)) for q in series if path("tol_" + q)}
    if not tols:
        raise ValueError("[post] necesita tol_cd o tol_cl")
    hit = num.first_settled(x, {q: series[q] for q in tols}, window, tols)
    j, lo = hit if hit else (n - 1, num.window_start(x, n - 1, window))

    out["historia"].update(unidad=unit, x_final=x[-1], t_final=t[-1])
    out["parada"] = {
        "criterio": "rango en %g %s: " % (window, unit)
                    + ", ".join("%s < %g" % (q, tol) for q, tol in tols.items()),
        "cumplio": "si" if hit else "no",
        "x_parada": x[j], "t_parada": t[j]}
    if not hit:
        rng = {q: max(series[q][lo:j + 1]) - min(series[q][lo:j + 1]) for q in tols}
        out["parada"]["nota"] = ("historial más corto que la ventana" if x[-1] - x[0] < window else
                                 "; ".join("%s rango %.3g" % (q, r) for q, r in rng.items()))
    for q, s in series.items():
        w = s[lo:j + 1]
        out["resultado"][q] = sum(w) / len(w)
        out["resultado"][q + "_rango"] = max(w) - min(w)
        out["resultado"][q + "_final"] = s[-1]
    return out


def cost(solver, logs, t_stop, procs, meta):
    tl, sec = parsers.timing(solver, logs[0]) if logs else ([], [])
    tb = meta.get("taubench_s") or None
    out = {"procesos": procs, "taubench_s": tb}
    if not sec:
        out["nota"] = "sin tiempos en el log"
        return out
    used = num.interp(t_stop, tl, sec) if t_stop is not None else sec[-1]
    out.update(segundos_solver=used, segundos_total=sec[-1], horas_nucleo=procs * used / 3600.0)
    if tb:
        out["costo_normalizado"] = procs * used / tb
    if t_stop is not None and t_stop > tl[-1]:
        out["nota"] = "el log termina antes que las fuerzas: tiempo recortado"
    return out


def write(rd, res):
    cfg.write_ini(os.path.join(rd, "result.ini"), res)
