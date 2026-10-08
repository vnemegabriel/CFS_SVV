"""Archivo de vuelo -> Cd en coasteo, contra el Cd(M) simulado.

Con acelerómetro (col_accel), que mide fuerza específica y no ve la gravedad:
    Cd = -m a / (0.5 rho V² A)                            a a lo largo del eje, alfa ~ 0
Sin acelerómetro, derivando la velocidad:
    Cd = -m (dV/dt + g sen(gamma)) / (0.5 rho V² A)      sen(gamma) = (dh/dt) / V

rho y c salen de la atmósfera estándar a la altura medida. V es la velocidad del
archivo: si es velocidad respecto del suelo, el viento entra como error.
"""
import csv
import datetime
import math
import os
import re
import zipfile

import num

G0, RGAS, GAMMA = 9.80665, 287.05, 1.4


def isa(h):
    """Densidad y velocidad del sonido, atmósfera estándar hasta 20 km."""
    if h < 11000.0:
        T = 288.15 - 0.0065 * h
        p = 101325.0 * (T / 288.15) ** 5.25588
    else:
        T = 216.65
        p = 22632.06 * math.exp(-G0 * (h - 11000.0) / (RGAS * T))
    return p / (RGAS * T), math.sqrt(GAMMA * RGAS * T)


def read_rows(path):
    if path.lower().endswith(".xlsx"):
        return _xlsx(path)
    with open(path, newline="", encoding="utf8", errors="replace") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        return list(csv.reader(f, dialect))


def _xlsx(path):
    z = zipfile.ZipFile(path)
    names = z.namelist()
    shared = []
    if "xl/sharedStrings.xml" in names:
        xml = z.read("xl/sharedStrings.xml").decode("utf8")
        shared = [re.sub(r"<[^>]+>", "", si) for si in re.findall(r"<si>(.*?)</si>", xml, re.S)]
    sheet = sorted(n for n in names if n.startswith("xl/worksheets/sheet"))[0]
    rows = []
    for body in re.findall(r"<row[^>]*>(.*?)</row>", z.read(sheet).decode("utf8"), re.S):
        row = {}
        for attrs, cell in re.findall(r"<c([^>]*?)(?:/>|>(.*?)</c>)", body, re.S):
            ref = re.search(r'r="([A-Z]+)\d+"', attrs)
            col = sum((ord(ch) - 64) * 26**i for i, ch in enumerate(reversed(ref.group(1)))) - 1
            v = re.search(r"<v>(.*?)</v>", cell or "")
            text = v.group(1) if v else re.sub(r"<[^>]+>", "", cell or "")
            row[col] = shared[int(text)] if 't="s"' in attrs and text else text
        rows.append([row.get(i, "") for i in range(max(row) + 1)] if row else [])
    return rows


def seconds(s):
    try:
        return float(s)
    except ValueError:
        return datetime.datetime.fromisoformat(s.strip()).timestamp()


def smooth(v, n):
    if n <= 1:
        return list(v)
    k = n // 2
    return [sum(v[max(0, i - k):i + k + 1]) / len(v[max(0, i - k):i + k + 1]) for i in range(len(v))]


def derivative(t, v):
    d = []
    for i in range(len(t)):
        a, b = max(0, i - 1), min(len(t) - 1, i + 1)
        d.append((v[b] - v[a]) / (t[b] - t[a]) if t[b] != t[a] else float("nan"))
    return d


def flight_cd(s, base_dir):
    """Devuelve (filas, nota). Filas vacías y nota explicativa si falta configuración."""
    need = ["file", "col_time", "col_alt", "col_speed", "mass_kg", "aref_m2"]
    missing = [k for k in need if not s.get(k)]
    if missing:
        return [], "telemetría: falta %s en [telemetry]" % ", ".join(missing)
    path = os.path.expanduser(s["file"])
    path = path if os.path.isabs(path) else os.path.join(base_dir, path)
    if not os.path.exists(path):
        return [], "telemetría: no existe %s" % path

    it, ih, iv = int(s["col_time"]), int(s["col_alt"]), int(s["col_speed"])
    ia = int(s["col_accel"]) if s.get("col_accel") else None
    scale = float(s.get("accel_scale") or 1.0)
    offset = float(s.get("alt_offset_m") or 0.0)
    fcol, fval = s.get("filter_col"), s.get("filter_value")
    data, skipped = [], 0
    for row in read_rows(path):
        try:
            if fcol and row[int(fcol)].strip() != fval:
                continue
            data.append((seconds(row[it]), float(row[ih]) + offset, float(row[iv]),
                         float(row[ia]) * scale if ia is not None else 0.0))
        except (IndexError, ValueError):
            skipped += 1
    if len(data) < 3:
        return [], "telemetría: menos de 3 filas legibles (%d descartadas)" % skipped

    t0 = data[0][0]
    lo, hi = float(s.get("t_start") or "-inf"), float(s.get("t_end") or "inf")
    data = [(t - t0, h, v, a) for t, h, v, a in data if lo <= t - t0 <= hi]
    t = [d[0] for d in data]
    n = int(s.get("smooth") or 1)
    h, v, acc = (smooth([d[k] for d in data], n) for k in (1, 2, 3))
    dv, dh = derivative(t, v), derivative(t, h)

    m, A = float(s["mass_kg"]), float(s["aref_m2"])
    rows = []
    for i in range(len(t)):
        if v[i] <= 1.0 or math.isnan(dv[i]):
            continue
        rho, c = isa(h[i])
        q = 0.5 * rho * v[i] ** 2 * A
        if ia is not None:
            cd = -m * acc[i] / q
        else:
            cd = -m * (dv[i] + G0 * max(-1.0, min(1.0, dh[i] / v[i]))) / q
        rows.append({"t": t[i], "h": h[i], "V": v[i], "mach": v[i] / c, "cd_vuelo": cd})
    method = "acelerómetro" if ia is not None else "derivada de la velocidad"
    if not rows:
        return [], "telemetría: ningún punto con V > 1 m/s"
    return rows, "telemetría (%s): %d puntos, Mach %.2f a %.2f, %d filas descartadas" % (
        method, len(rows), min(r["mach"] for r in rows), max(r["mach"] for r in rows), skipped)


def compare(rows, machs, cds):
    for r in rows:
        sim = num.pchip(machs, cds, r["mach"]) if len(machs) >= 2 else None
        r["cd_sim"] = sim
        r["diferencia"] = None if sim is None else r["cd_vuelo"] - sim
    return rows
