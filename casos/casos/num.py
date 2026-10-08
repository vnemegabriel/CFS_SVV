"""Numérica sin dependencias: interpolación, criterio de parada, PCHIP, Richardson."""
import bisect
import math


def interp(x, xs, ys):
    """Lineal, recortada a los extremos."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    i = bisect.bisect_right(xs, x)
    w = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
    return ys[i - 1] + w * (ys[i] - ys[i - 1])


def window_start(x, j, window):
    return bisect.bisect_left(x, x[j] - window, 0, j + 1)


def settled(values, tol):
    if any(math.isnan(v) or math.isinf(v) for v in values):
        return False
    return max(values) - min(values) < tol


def first_settled(x, series, window, tols):
    """Primer índice j donde cada serie varía menos que su tolerancia en [x_j - window, x_j].

    Devuelve (j, inicio de la ventana) o None.
    """
    for j in range(len(x)):
        if x[j] - x[0] < window:
            continue
        lo = window_start(x, j, window)
        if all(settled(s[lo:j + 1], tols[q]) for q, s in series.items()):
            return j, lo
    return None


def pchip(xs, ys, xq):
    """Hermite cúbica monótona (Fritsch–Carlson). None fuera del rango de datos."""
    n = len(xs)
    if n == 0 or xq < xs[0] or xq > xs[-1]:
        return None
    if n == 1:
        return ys[0]
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    d = [(ys[i + 1] - ys[i]) / h[i] for i in range(n - 1)]
    m = [0.0] * n
    if n == 2:
        m = [d[0], d[0]]
    else:
        for i in range(1, n - 1):
            if d[i - 1] * d[i] > 0:
                w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
                m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
        m[0] = _end_slope(h[0], h[1], d[0], d[1])
        m[-1] = _end_slope(h[-1], h[-2], d[-1], d[-2])
    i = min(max(bisect.bisect_right(xs, xq) - 1, 0), n - 2)
    t = (xq - xs[i]) / h[i]
    return ((2 * t**3 - 3 * t**2 + 1) * ys[i] + (t**3 - 2 * t**2 + t) * h[i] * m[i]
            + (-2 * t**3 + 3 * t**2) * ys[i + 1] + (t**3 - t**2) * h[i] * m[i + 1])


def _end_slope(h0, h1, d0, d1):
    m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if m * d0 <= 0:
        return 0.0
    if d0 * d1 < 0 and abs(m) > abs(3 * d0):
        return 3 * d0
    return m


def richardson(q3, q2, q1, r=2.0):
    """q3 gruesa, q1 fina. Devuelve (dict p/qinf/u, nota). dict es None si no aplica."""
    e32, e21 = q3 - q2, q2 - q1
    if e21 == 0 or e32 == 0:
        return None, "sin variación entre niveles"
    if e32 * e21 < 0:
        return None, "convergencia no monótona"
    p = math.log(e32 / e21) / math.log(r)
    if p <= 0:
        return None, "orden observado %.2f: no converge" % p
    k = r**p - 1
    return {"p": p, "qinf": q1 + (q1 - q2) / k, "u": 1.25 * abs(q1 - q2) / k}, ""
