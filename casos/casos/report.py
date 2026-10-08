"""Tablas y figuras de un caso. Se rehacen enteras desde los archivos de cada corrida."""
import collections
import csv
import math
import os
import time

import cfg
import num
import runner
import telemetry

ORDER = ("id", "estado", "motivo", "raiz", "matriz.", "corrida.", "malla.", "historia.",
         "parada.", "resultado.", "errores.", "costo.")
KINDS = ("cost_curve", "cd_mach", "errores")


def key_order(k):
    for i, prefix in enumerate(ORDER):
        if k == prefix or (prefix.endswith(".") and k.startswith(prefix)):
            return i, k
    return len(ORDER), k


def write_csv(path, rows, cols=None):
    cols = cols or sorted({k for r in rows for k in r}, key=key_order)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, cols, extrasaction="ignore")
        w.writeheader()
        w.writerows({k: cfg.fmt(v) for k, v in r.items()} for r in rows)
    os.replace(tmp, path)


def col(row, key):
    for k in ("matriz." + key, "corrida." + key, key):
        if row.get(k, "") != "":
            return row[k]
    return ""


def gather(caso, roots, stale_s=None):
    stale_s = stale_s or 5 * runner.limits(caso)[2]
    rows = []
    for m in caso.rows():
        base = {"id": m["id"], "caso": caso.name}
        base.update({"matriz." + k: v for k, v in m.items() if k != "id"})
        dirs = [os.path.join(r, caso.name, m["id"]) for r in roots]
        dirs = [d for d in dirs if os.path.exists(os.path.join(d, "status"))]
        if not dirs:
            rows.append(dict(base, estado="pendiente", motivo="", raiz=""))
        for rd in dirs:
            st = runner.read_status(rd)
            row = dict(base, estado=st["estado"], motivo=st["motivo"],
                       raiz=os.path.dirname(os.path.dirname(rd)))
            age = runner.heartbeat_age(rd)
            if st["estado"] == "corriendo" and (age is None or age > stale_s):
                row["estado"] = "corriendo?"
                row["motivo"] = "sin latido del runner; ¿se cortó? (`casos.py unlock`)"
            row.update(cfg.read_ini_flat(os.path.join(rd, "result.ini")))
            rows.append(row)
    return rows


def regenerate(caso, roots):
    out = caso.results_dir
    os.makedirs(out, exist_ok=True)
    notes = []
    rows = gather(caso, roots)
    for name in (caso.get("report", "join") or "").split():
        try:
            rows += gather(cfg.Caso(name), roots)
        except Exception as e:
            notes.append("join %s: %s" % (name, e))
    write_csv(os.path.join(out, "results.csv"), rows)

    kind = caso.get("report", "kind")
    try:
        if kind == "cost_curve":
            cost_curve(caso, rows, out, notes)
        elif kind == "cd_mach":
            cd_mach(caso, rows, out, notes)
        elif kind == "errores":
            error_table(caso, rows, out, notes)
        elif kind:
            notes.append("[report] kind '%s' desconocido (%s)" % (kind, ", ".join(KINDS)))
    except Exception as e:
        notes.append("reporte %s falló: %s: %s" % (kind, type(e).__name__, e))

    no_tb = sorted({r.get("corrida.maquina", "") for r in rows
                    if r.get("costo.segundos_solver") and not r.get("costo.taubench_s")})
    if no_tb:
        notes.append("costo normalizado vacío en %s: falta taubench_s en maquinas.ini" % ", ".join(no_tb))

    count = collections.Counter(r["estado"] for r in rows)
    text = ["caso %s   %s" % (caso.name, time.strftime("%Y-%m-%d %H:%M")),
            "raíces: " + " ".join(roots),
            "   ".join("%s %d" % kv for kv in sorted(count.items())), ""]
    bad = [r for r in rows if r["estado"] not in ("ok", "pendiente")]
    if bad:
        text.append("para revisar:")
        text += ["  %-34s %-14s %s" % (r["id"], r["estado"], r["motivo"]) for r in bad]
        text.append("")
    if notes:
        text.append("notas:")
        text += ["  - " + n for n in notes]
    text = "\n".join(text) + "\n"
    with open(os.path.join(out, "resumen.txt"), "w", encoding="utf8") as f:
        f.write(text)
    return text


# ---------------------------------------------------------------- curva costo–precisión

def cost_curve(caso, rows, out, notes):
    s = dict(caso.items("report"))
    group = s.get("group", "mach").split()
    family = s.get("family", "solver").split()
    level, order = s.get("level", "level"), s.get("order", "P")
    q, ratio = s.get("quantity", "cd"), float(s.get("ratio", "2"))
    reference = dict(kv.split("=", 1) for kv in s.get("reference", "").split())
    if not reference:
        notes.append("[report] sin reference: la curva queda sin error")

    done = [r for r in rows if r["estado"] == "ok" and r.get("resultado." + q, "") != ""]
    gkey = lambda r: tuple(col(r, k) for k in group)

    refs, ref_rows = {}, []
    for g in sorted({gkey(r) for r in done}):
        line = dict(zip(group, g))
        levels = collections.defaultdict(dict)
        for r in done:
            if gkey(r) == g and reference and all(col(r, k) == v for k, v in reference.items()):
                levels[col(r, order)][int(col(r, level))] = float(r["resultado." + q])
        usable = {p: lv for p, lv in levels.items() if len(lv) >= 3}
        if not reference:
            continue
        if not usable:
            line["nota"] = "menos de 3 niveles ok en la familia de referencia"
        else:
            p = max(usable, key=float)
            finest = sorted(usable[p])[-3:]
            rich, note = num.richardson(*(usable[p][k] for k in finest), r=ratio)
            line.update(orden=p, niveles=" ".join(map(str, finest)), nota=note)
            if rich:
                line.update(orden_observado=rich["p"], referencia=rich["qinf"], incertidumbre=rich["u"])
                refs[g] = rich
        ref_rows.append(line)
    if ref_rows:
        write_csv(os.path.join(out, "referencia.csv"), ref_rows,
                  group + ["orden", "niveles", "orden_observado", "referencia", "incertidumbre", "nota"])

    curve = []
    for r in done:
        ref, value = refs.get(gkey(r)), float(r["resultado." + q])
        curve.append(dict(zip(group, gkey(r)),
                          familia=" ".join(col(r, k) for k in family),
                          nivel=col(r, level), orden=col(r, order),
                          dof=r.get("malla.dof", ""), h_eff=r.get("malla.h_eff", ""),
                          horas_nucleo=r.get("costo.horas_nucleo", ""),
                          costo_normalizado=r.get("costo.costo_normalizado", ""),
                          valor=value,
                          error=abs(value - ref["qinf"]) if ref else "",
                          incertidumbre_ref=ref["u"] if ref else ""))
    write_csv(os.path.join(out, "curva.csv"), curve,
              group + ["familia", "nivel", "orden", "dof", "h_eff", "horas_nucleo",
                       "costo_normalizado", "valor", "error", "incertidumbre_ref"])
    figures_cost(out, curve, group, q, notes)


# ---------------------------------------------------------------- Cd contra Mach

def cd_mach(caso, rows, out, notes):
    s = dict(caso.items("report"))
    family, mcol = s.get("family", "solver").split(), s.get("mach", "mach")
    done = [r for r in rows if r["estado"] == "ok" and r.get("resultado.cd", "") != ""]
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in done:
        by[" ".join(col(r, k) for k in family)][float(col(r, mcol))].append(float(r["resultado.cd"]))

    points, dense, curves = [], [], {}
    for fam, pts in sorted(by.items()):
        ms = sorted(pts)
        cds = [sum(pts[m]) / len(pts[m]) for m in ms]
        curves[fam] = (ms, cds)
        points += [dict(familia=fam, mach=m, cd=c, corridas=len(pts[m])) for m, c in zip(ms, cds)]
        if len(ms) >= 2:
            n = int(round((ms[-1] - ms[0]) / 0.01))
            dense += [dict(familia=fam, mach=ms[0] + i * (ms[-1] - ms[0]) / n,
                           cd=num.pchip(ms, cds, ms[0] + i * (ms[-1] - ms[0]) / n)) for i in range(n + 1)]
    write_csv(os.path.join(out, "cd_mach.csv"), points, ["familia", "mach", "cd", "corridas"])
    write_csv(os.path.join(out, "cd_mach_interp.csv"), dense, ["familia", "mach", "cd"])
    if len(points) < 5:
        notes.append("Cd(M): %d puntos ok; el objetivo son al menos 5 entre Mach 0 y 1.8" % len(points))

    flight = []
    if caso.ini.has_section("telemetry"):
        tele = dict(caso.items("telemetry"))
        flight, note = telemetry.flight_cd(tele, caso.dir)
        notes.append(note)
        fam = tele.get("family") or (sorted(curves)[0] if curves else "")
        ms, cds = curves.get(fam, ([], []))
        if flight:
            if len(ms) < 2:
                notes.append("telemetría: '%s' tiene menos de 2 Mach ok, sin comparación" % fam)
            telemetry.compare(flight, ms, cds)
            write_csv(os.path.join(out, "telemetria_cd.csv"), flight,
                      ["t", "h", "V", "mach", "cd_vuelo", "cd_sim", "diferencia"])
    figures_cd(out, curves, flight, notes)


# ---------------------------------------------------------------- normas de error y criterios

def norms(row, name):
    prefix = "errores.%s_" % name
    return {k[len(prefix):]: float(v) for k, v in row.items() if k.startswith(prefix) and v != ""}


def criterion(caso, row, done, columns):
    """Evalúa [criteria] para el caso de la fila. Devuelve (si|no|pendiente|"", detalle).

    linf <= X                  todas las normas L∞ bajo X
    igual stab=dg tol=X        normas L2 iguales a las de la fila gemela con stab=dg
    orden h variable=V margen=X   orden observado en h de la norma L2 de V >= P+1-X
    """
    spec = caso.get("criteria", col(row, caso.get("report", "case", "caso")))
    if not spec:
        return "", ""
    words = spec.split()
    opts = dict(w.split("=", 1) for w in words[1:] if "=" in w)

    if words[0] == "linf":
        values = norms(row, "Linf")
        if not values:
            return "pendiente", "sin normas L∞"
        worst = max(values.values())
        return ("si" if worst <= float(words[-1]) else "no"), "L∞ máx %.3g" % worst

    if words[0] == "igual":
        ref = {k: v for k, v in opts.items() if k != "tol"}
        if all(col(row, k) == v for k, v in ref.items()):
            return "", "es la referencia"
        twins = [t for t in done if all(col(t, c) == ref.get(c, col(row, c)) for c in columns)]
        if not twins:
            return "pendiente", "falta la corrida con " + " ".join("%s=%s" % kv for kv in ref.items())
        a, b = norms(row, "L2"), norms(twins[0], "L2")
        diff = max(abs(a[k] - b[k]) / max(abs(b[k]), 1e-300) for k in b if k in a)
        return ("si" if diff <= float(opts.get("tol", "1e-10")) else "no"), "diferencia relativa %.3g" % diff

    if words[0] == "orden":
        hcol, var = words[1], opts.get("variable", "rho")
        same = [t for t in done if all(col(t, c) == col(row, c) for c in columns if c != hcol)]
        pts = sorted({(float(col(t, hcol)), norms(t, "L2").get(var)) for t in same
                      if col(t, hcol) != "" and norms(t, "L2").get(var)}, reverse=True)
        if len(pts) < 2:
            return "pendiente", "hacen falta al menos dos valores de %s" % hcol
        (h1, e1), (h2, e2) = pts[-2], pts[-1]
        p = math.log(e1 / e2) / math.log(h1 / h2)
        formal = float(col(row, "P")) + 1
        return ("si" if p >= formal - float(opts.get("margen", "0.3")) else "no"), \
            "orden %.2f, formal %g" % (p, formal)

    return "pendiente", "criterio desconocido: " + spec


def error_table(caso, rows, out, notes):
    columns = caso.get("report", "columns", "caso stab P").split()
    done = [r for r in rows if r["estado"] == "ok"]
    names = sorted({k[len("errores."):] for r in done for k in r if k.startswith("errores.")})
    table = []
    for r in done:
        line = {c: col(r, c) for c in columns}
        line.update({n: r["errores." + n] for n in names if r.get("errores." + n, "") != ""})
        line["cumple"], line["detalle"] = criterion(caso, r, done, columns)
        table.append(line)
    write_csv(os.path.join(out, "errores.csv"), table, columns + names + ["cumple", "detalle"])
    count = collections.Counter(t["cumple"] for t in table if t["cumple"])
    if count:
        notes.append("criterios: " + "   ".join("%s %d" % kv for kv in sorted(count.items())))
    for t in table:
        if t["cumple"] == "no":
            notes.append("NO CUMPLE %s: %s" % (" ".join(t[c] for c in columns), t["detalle"]))


# ---------------------------------------------------------------- figuras

def pyplot(notes):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        notes.append("sin matplotlib: no se generan figuras")
        return None


def figures_cost(out, curve, group, q, notes):
    plt = pyplot(notes)
    rows = [r for r in curve if r["error"] != "" and r["error"] > 0]
    if not plt or not rows:
        return
    for g in sorted({tuple(r[k] for k in group) for r in rows}):
        sel = [r for r in rows if tuple(r[k] for k in group) == g]
        x = "costo_normalizado" if all(r["costo_normalizado"] != "" for r in sel) else "horas_nucleo"
        fig, ax = plt.subplots(figsize=(6, 4.5))
        for fam in sorted({r["familia"] for r in sel}):
            pts = sorted((float(r[x]), r["error"]) for r in sel if r["familia"] == fam and r[x] != "")
            if pts:
                ax.loglog(*zip(*pts), "o-", label=fam)
        ax.set_xlabel(x.replace("_", " "))
        ax.set_ylabel("|%s − referencia|" % q)
        ax.set_title(" ".join("%s=%s" % kv for kv in zip(group, g)))
        ax.grid(alpha=0.3, which="both")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, "curva_%s.png" % "_".join(g)), dpi=130)
        plt.close(fig)


def figures_cd(out, curves, flight, notes):
    plt = pyplot(notes)
    if not plt or not curves:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for fam, (ms, cds) in sorted(curves.items()):
        line, = ax.plot(ms, cds, "o", label=fam)
        if len(ms) >= 2:
            xs = [ms[0] + i * (ms[-1] - ms[0]) / 200 for i in range(201)]
            ax.plot(xs, [num.pchip(ms, cds, x) for x in xs], "-", color=line.get_color())
    if flight:
        ax.plot([r["mach"] for r in flight], [r["cd_vuelo"] for r in flight], ".", color="0.5",
                ms=3, label="vuelo")
    ax.set_xlabel("Mach")
    ax.set_ylabel("Cd")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "cd_mach.png"), dpi=130)
    plt.close(fig)
