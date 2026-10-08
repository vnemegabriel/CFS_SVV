#!/usr/bin/env python3
"""Validación regularizada (RRSI, Xia et al. 2026, arXiv 2609.24972) para la calibración de SVV.

Un candidato = commit de ~/nektar (rama feature/svv-cfs) + parámetros de sesión (-P clave=valor).
Controles de RRSI adaptados:
  propuesta:  presupuesto de ediciones con calendario coseno; historial de intentos; poda de componentes.
  selección:  filtro de filtración (regex + crítico opcional); tolerancia de ruido medida; piso
              mejor-τ por caso; regla costo–precisión; holdout separado que el proponente no ve.

Subcomandos:
  ruido         corre el incumbente N veces y fija τ
  ronda         evalúa un candidato sobre 'evolve' y decide
  holdout       evalúa el incumbente sobre 'holdout' (sólo orquestador)
  transparencia dev vs base sin SVV: normas L2 idénticas
  historial     últimas entradas del ledger, una línea cada una
  podar         componentes sin ganancia en la ventana de poda
  presupuesto   ediciones permitidas en la próxima ronda
"""
import argparse, concurrent.futures as cf, json, math, os, re, shutil, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
H = os.path.dirname(HERE)
NEK = os.path.expanduser("~/nektar")
CONJ = os.path.join(H, "validacion", "conjuntos.json")
LEDGER = os.path.join(H, "validacion", "ledger.jsonl")
HOLD = os.path.join(H, "validacion", "holdout.jsonl")
STATE = os.path.join(H, "validacion", "estado.json")
RUNS = os.path.expanduser("~/runs-dev/validar")

RE_L2 = re.compile(r"L 2 error \(variable (\w+)\)\s*:\s*([-+0-9.eE]+|nan|inf)", re.I)
RE_T = re.compile(r"Total Computation Time\s*[=:]\s*([0-9.eE+-]+)")


def cargar(p, defecto):
    try:
        with open(p) as f:
            return json.load(f)
    except FileNotFoundError:
        return defecto


def guardar(p, obj):
    with open(p, "w") as f:
        json.dump(obj, f, indent=1)


def ledger(p=LEDGER):
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return [json.loads(l) for l in f if l.strip()]


def anexar(p, e):
    with open(p, "a") as f:
        f.write(json.dumps(e, ensure_ascii=False) + "\n")


def git(*a):
    return subprocess.run(["git", "-C", NEK, *a], capture_output=True, text=True).stdout.strip()


# ---------- corrida de un caso ----------

def correr(caso, inst, params, tag):
    """Devuelve {'score','costo','l2','ok'}; score = -log10(L2 de la variable del caso)."""
    d = os.path.join(RUNS, tag, caso["id"])
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    archivos = [os.path.expanduser(a) for a in caso["sesion"]]
    for a in archivos:
        shutil.copy(a, d)
    pv = dict(caso.get("params", {}))
    pv.update(params)
    cmd = [os.path.join(HERE, "nk"), inst]
    if caso.get("procs", 1) > 1:
        cmd += ["mpirun", "--bind-to", "core", "-np", str(caso["procs"])]
    cmd += ["CompressibleFlowSolver", *[os.path.basename(a) for a in archivos]]
    for k, v in pv.items():
        cmd += ["-P", f"{k}={v}"]
    env = dict(os.environ, OMPI_MCA_hwloc_base_use_hwthreads_as_cpus="1")
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env,
                           timeout=caso.get("timeout_s", 1800))
        out, rc = r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired:
        out, rc = "TIMEOUT", -1
    with open(os.path.join(d, "log"), "w") as f:
        f.write(" ".join(cmd) + "\n" + out)
    l2 = {m.group(1): float(m.group(2)) for m in RE_L2.finditer(out)}
    mt = RE_T.search(out)
    costo = float(mt.group(1)) if mt else time.time() - t0
    v = l2.get(caso.get("var", "rho"))
    ok = rc == 0 and v is not None and math.isfinite(v) and v > 0
    return {"ok": ok, "score": -math.log10(v) if ok else None, "costo": costo, "l2": l2,
            "rc": rc, "log": os.path.join(d, "log")}


def correr_conjunto(nombre, inst, params, tag, paralelo):
    casos = cargar(CONJ, {})[nombre]
    faltan = [f"{c['id']}: {a}" for c in casos for a in c["sesion"] if not os.path.exists(os.path.expanduser(a))]
    pend = [f"{c['id']}: {k}" for c in casos for k, v in {**c.get("params", {}), **params}.items() if v == "PENDIENTE"]
    if faltan or pend:
        sys.exit(f"{nombre}: no se puede correr.\n  faltan sesiones: {faltan}\n  PENDIENTE: {pend}")
    if paralelo:
        with cf.ThreadPoolExecutor(max_workers=paralelo) as ex:
            res = list(ex.map(lambda c: correr(c, inst, params, tag), casos))
    else:
        res = [correr(c, inst, params, tag) for c in casos]
    return {c["id"]: r for c, r in zip(casos, res)}


# ---------- controles RRSI ----------

def presupuesto(ronda, cfg):
    kmax, kmin, T = cfg["k_max"], cfg["k_min"], cfg["rondas"]
    t = min(ronda, T)
    return int(round(kmin + 0.5 * (kmax - kmin) * (1 + math.cos(math.pi * t / T))))


def n_ediciones(est, params):
    base = est.get("commit")
    archivos = set(git("diff", "--name-only", base).splitlines()) if base else set()
    archivos |= set(l[3:] for l in git("status", "--porcelain").splitlines())
    pv = est.get("params", {})
    cambiados = {k for k in set(pv) | set(params) if str(pv.get(k)) != str(params.get(k))}
    return len(archivos), cambiados


def filtracion(est, cfg):
    base = est.get("commit") or git("merge-base", "HEAD", "master")
    diff = git("diff", base, "--", "solvers/CompressibleFlowSolver")
    agregadas = [l[1:] for l in diff.splitlines() if l.startswith("+") and not l.startswith("+++")]
    codigo = [l for l in agregadas if not l.strip().startswith("//")]
    hallazgos = [l.strip() for l in codigo if re.search(cfg["regex_filtracion"], l, re.I)]
    neto = len(agregadas) - sum(1 for l in diff.splitlines() if l.startswith("-") and not l.startswith("---"))
    return hallazgos, neto, diff


def critico(diff, hipotesis):
    tarea = (f"Hipótesis declarada: {hipotesis}\n\nDiff (git -C ~/nektar diff):\n```diff\n{diff[:60000]}\n```")
    r = subprocess.run([os.path.join(HERE, "oc"), "svv-critic", tarea], capture_output=True, text=True)
    m = re.findall(r"\{.*\"verdict\".*\}", r.stdout)
    try:
        return json.loads(m[-1])
    except Exception:
        return {"verdict": "NO_ANSWER", "reasons": [r.stderr[-300:]]}


def decidir(res, est, ruido, cfg, neto):
    tau, tau_c = ruido["tau"], ruido["tau_costo_rel"]
    inc, mejor = est["scores"], est["mejor"]
    fallos = [k for k, r in res.items() if not r["ok"]]
    if fallos:
        return "REJECT", f"corridas fallidas: {fallos}", None, None
    bajo_piso = [k for k, r in res.items() if r["score"] < mejor.get(k, -1e9) - tau]
    if bajo_piso:
        return "REJECT", f"bajo el piso mejor-τ en {bajo_piso}", None, None
    d = statistics.mean(r["score"] - inc[k] for k, r in res.items())
    c_new = sum(r["costo"] for r in res.values())
    c_old = sum(est["costos"].values())
    dc = (c_new - c_old) / c_old
    if d > tau:
        if dc <= cfg["lambda_costo"] * d + tau_c:
            return "ACCEPT", f"Δ={d:+.3f} > τ, costo {dc:+.1%} dentro de λ·Δ", d, dc
        return "REJECT", f"Δ={d:+.3f} pero costo {dc:+.1%} > λ·Δ", d, dc
    if d >= -tau and (dc < -tau_c or neto < 0):
        return "ACCEPT", f"Δ={d:+.3f} en banda de ruido con simplificación (costo {dc:+.1%}, neto {neto} líneas)", d, dc
    return "REJECT", f"Δ={d:+.3f} sin evidencia sobre τ={tau:.3g}", d, dc


# ---------- subcomandos ----------

def cmd_ruido(a):
    cfg = cargar(CONJ, {})["rrsi"]
    est = cargar(STATE, {"params": {}})
    est.setdefault("params", {})
    est["commit"] = est.get("commit") or git("rev-parse", "HEAD")
    reps = [correr_conjunto("evolve", "dev", est["params"], f"ruido-{i}", a.paralelo)
            for i in range(a.n or cfg["repeticiones_ruido"])]
    ids = reps[0].keys()
    if not all(r[k]["ok"] for r in reps for k in ids):
        sys.exit("ruido: alguna corrida del incumbente falló; ver logs en ~/runs-dev/validar/ruido-*")
    tau = max(max(r[k]["score"] for r in reps) - min(r[k]["score"] for r in reps) for k in ids)
    tot = [sum(r[k]["costo"] for k in ids) for r in reps]
    ruido = {"tau": max(tau, cfg["tau_min"]),
             "tau_costo_rel": statistics.pstdev(tot) / statistics.mean(tot) if len(tot) > 1 else 0.05}
    est.update(scores={k: statistics.mean(r[k]["score"] for r in reps) for k in ids},
               costos={k: statistics.mean(r[k]["costo"] for r in reps) for k in ids}, ruido=ruido)
    est.setdefault("mejor", dict(est["scores"]))
    est.setdefault("ronda", 0)
    guardar(STATE, est)
    print(f"τ={ruido['tau']:.3g}  τ_costo={ruido['tau_costo_rel']:.1%}  incumbente {est['commit'][:9]} {est['params']}")


def cmd_ronda(a):
    cj = cargar(CONJ, {})
    cfg = cj["rrsi"]
    est = cargar(STATE, None)
    if not est or "ruido" not in est:
        sys.exit("ronda: correr primero `validar.py ruido` para fijar incumbente y τ")
    params = dict(kv.split("=", 1) for kv in a.params)
    ronda = est["ronda"] + 1
    k = presupuesto(ronda, cfg)
    n_arch, n_par = n_ediciones(est, params)
    ent = {"ronda": ronda, "fecha": time.strftime("%F %T"), "commit": git("rev-parse", "HEAD"),
           "params": params, "hipotesis": a.hipotesis, "componentes": a.componentes,
           "ediciones": n_arch + len(n_par), "presupuesto": k}
    def fin(veredicto, motivo, **extra):
        ent.update(veredicto=veredicto, motivo=motivo, **extra)
        est["ronda"] = ronda
        guardar(STATE, est)
        anexar(LEDGER, ent)
        print(f"ronda {ronda}: {veredicto} — {motivo}")
        sys.exit(0 if veredicto == "ACCEPT" else 1)

    if n_arch + len(n_par) > k:
        fin("REJECT", f"{n_arch} archivos + {len(n_par)} parámetros > presupuesto {k}")
    hall, neto, diff = filtracion(est, cfg)
    if hall:
        fin("LEAK", "regex: " + " | ".join(hall[:3]))
    if a.critico and diff:
        v = critico(diff, a.hipotesis)
        ent["critico"] = v
        if v.get("verdict") in ("LEAK", "INERT"):
            fin(v["verdict"], "crítico: " + "; ".join(v.get("reasons", [])[:2]))
    res = correr_conjunto("evolve", "dev", params, f"ronda-{ronda}", a.paralelo)
    ent["scores"] = {k_: r["score"] for k_, r in res.items()}
    ent["costos"] = {k_: r["costo"] for k_, r in res.items()}
    v, motivo, d, dc = decidir(res, est, est["ruido"], cfg, neto)
    ent.update(delta=d, delta_costo=dc)
    if v == "ACCEPT":
        est.update(commit=ent["commit"], params=params, scores=ent["scores"], costos=ent["costos"])
        est["mejor"] = {k_: max(est["mejor"].get(k_, -1e9), s) for k_, s in ent["scores"].items()}
    fin(v, motivo)


def cmd_holdout(a):
    est = cargar(STATE, None) or sys.exit("holdout: no hay incumbente")
    res = correr_conjunto("holdout", "dev", est.get("params", {}), "holdout", a.paralelo)
    base = correr_conjunto("holdout", "base", {}, "holdout-base", a.paralelo) if a.base else {}
    e = {"fecha": time.strftime("%F %T"), "commit": est["commit"], "params": est.get("params"),
         "scores": {k: r["score"] for k, r in res.items()},
         "costos": {k: r["costo"] for k, r in res.items()},
         "base": {k: r["score"] for k, r in base.items()}}
    anexar(HOLD, e)
    for k, r in res.items():
        b = f"  base {base[k]['score']:.3f}" if k in base and base[k]["ok"] else ""
        print(f"{k:28s} {'%.3f' % r['score'] if r['ok'] else 'FALLA'}  {r['costo']:.1f}s{b}")


def cmd_transparencia(a):
    if not os.path.exists(os.path.join(NEK, "build-svv/dist/bin/CompressibleFlowSolver")):
        sys.exit("transparencia: falta el build dev (harness/bin/svv-build)")
    cj = cargar(CONJ, {})
    casos = cj["transparencia"]
    peor = 0.0
    for c in casos:
        rb, rd = correr(c, "base", {}, "transp-base"), correr(c, "dev", {}, "transp-dev")
        for var, eb in rb["l2"].items():
            ed = rd["l2"].get(var)
            rel = abs(ed - eb) / max(abs(eb), 1e-300) if ed is not None else float("inf")
            peor = max(peor, rel)
            if rel > cj["rrsi"]["tol_transparencia"]:
                print(f"{c['id']} {var}: base {eb:.10e} dev {ed}")
    print(f"transparencia: {'OK' if peor <= cj['rrsi']['tol_transparencia'] else 'FALLA'} (máx rel {peor:.2e})")
    sys.exit(0 if peor <= cj["rrsi"]["tol_transparencia"] else 1)


def cmd_historial(a):
    for e in ledger()[-a.n:]:
        d = f"Δ={e['delta']:+.3f}" if e.get("delta") is not None else ""
        print(f"r{e['ronda']} {e['veredicto']} [{','.join(e.get('componentes') or [])}] {d} — {e.get('hipotesis','')[:70]} ({e['motivo'][:60]})")


def cmd_podar(a):
    cfg = cargar(CONJ, {})["rrsi"]
    L = ledger()[-cfg["ventana_poda"]:]
    vistos, utiles = set(), set()
    for e in L:
        for c in e.get("componentes") or []:
            vistos.add(c)
            if e["veredicto"] == "ACCEPT" and (e.get("delta") or 0) > 0:
                utiles.add(c)
    estancado = len(L) >= cfg["estancamiento"] and all(
        e["veredicto"] != "ACCEPT" or abs(e.get("delta") or 0) <= cargar(STATE, {}).get("ruido", {}).get("tau", 0)
        for e in L[-cfg["estancamiento"]:])
    print("podar:", ", ".join(sorted(vistos - utiles)) or "nada")
    if estancado:
        print("estancado: reservar la próxima propuesta para un componente no explorado")


def cmd_presupuesto(a):
    cfg = cargar(CONJ, {})["rrsi"]
    est = cargar(STATE, {"ronda": 0})
    print(presupuesto(est.get("ronda", 0) + 1, cfg))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    for n in ("ruido", "ronda", "holdout"):
        q = s.add_parser(n)
        q.add_argument("--paralelo", type=int, default=0, help="casos simultáneos (afecta la medición de costo)")
    s.choices["ruido"].add_argument("-n", type=int)
    r = s.choices["ronda"]
    r.add_argument("--hipotesis", required=True)
    r.add_argument("--componentes", type=lambda x: x.split(","), default=[])
    r.add_argument("--params", nargs="*", default=[], metavar="CLAVE=VALOR")
    r.add_argument("--critico", action="store_true", help="segunda opinión de oc svv-critic")
    s.choices["holdout"].add_argument("--base", action="store_true", help="también corre base")
    s.add_parser("transparencia")
    s.add_parser("historial").add_argument("-n", type=int, default=10)
    s.add_parser("podar")
    s.add_parser("presupuesto")
    a = p.parse_args()
    globals()["cmd_" + a.cmd](a)


if __name__ == "__main__":
    main()
