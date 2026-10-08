#!/usr/bin/env python3
"""Read a Nektar++ session: PARAMETERS, NUMMODES, and mesh metrics.

Everything downstream (dt estimate, force coefficients, FieldConvert expressions)
derives from here, so changing Mach / Reynolds / alpha in conditions.xml is the
only edit needed -- nothing is hard-coded in the post-processing.
"""
import math, re, sys

# Nektar++ expression functions available inside <P> ... </P>
_NS = {
    "sqrt": math.sqrt, "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "asin": math.asin, "acos": math.acos, "atan": math.atan, "atan2": math.atan2,
    "exp": math.exp, "log": math.log, "abs": abs, "pow": pow,
    "PI": math.pi, "ceil": math.ceil, "floor": math.floor,
    "max": max, "min": min,
}


def parameters(path="conditions.xml"):
    """Evaluate every <P> name = expr </P> in declaration order."""
    txt = open(path, encoding="utf8", errors="ignore").read()
    txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    env = dict(_NS)
    out = {}
    for name, expr in re.findall(r"<P>\s*([A-Za-z_]\w*)\s*=\s*([^<]+?)\s*</P>", txt):
        try:
            v = eval(expr, {"__builtins__": {}}, env)
        except Exception:
            continue                      # non-numeric / forward reference
        env[name] = out[name] = v
    return out


def solverinfo(path="conditions.xml"):
    """Every <I PROPERTY="..." VALUE="..."/> in the session, as a dict."""
    txt = open(path, encoding="utf8", errors="ignore").read()
    txt = re.sub(r"<!--.*?-->", "", txt, flags=re.S)
    return {k: v for k, v in re.findall(
        r'<I\s+PROPERTY\s*=\s*"([^"]+)"\s+VALUE\s*=\s*"([^"]+)"', txt)}


def tag(path="conditions.xml"):
    """Short label for output files: 'euler' or 'ns'."""
    eq = solverinfo(path).get("EQType", "")
    return "euler" if eq.startswith("Euler") else "ns"


def nummodes(path="conditions.xml"):
    txt = open(path, encoding="utf8", errors="ignore").read()
    m = re.search(r'NUMMODES\s*=\s*"(\d+)', txt)
    return int(m.group(1)) if m else None


def mesh_metrics(path):
    """Minimum straight edge length, element count, and wall/farfield edge counts."""
    s = open(path, encoding="utf8", errors="ignore").read()
    V = {int(v): tuple(map(float, c.split()[:2]))
         for v, c in re.findall(r'<V ID="(\d+)">\s*([-\d.eE+ \t]+?)</V>', s)}
    E = {int(e): tuple(map(int, c.split()))
         for e, c in re.findall(r'<E ID="(\d+)">\s*(\d+\s+\d+)\s*</E>', s)}
    L = [math.dist(V[a], V[b]) for a, b in E.values()]
    nel = len(re.findall(r"<[QTHPRA] ID=", s))
    return {"hmin": min(L), "hmax": max(L), "nedge": len(L), "nelmt": nel}


def reference(params):
    """Freestream reference quantities used for non-dimensionalisation."""
    g, R = params["Gamma"], params.get("GasConstant", 287.058)
    rho, p = params["rhoInf"], params["pInf"]
    T = p / (R * rho)
    c = math.sqrt(g * R * T)
    M = params["Mach"]
    U = M * c
    L = params.get("Lref", 1.0)
    mu = params["mu"]
    return {"gamma": g, "R": R, "rho": rho, "p": p, "T": T, "c": c, "M": M,
            "U": U, "L": L, "mu": mu, "nu": mu / rho, "Re": rho * U * L / mu,
            "q": 0.5 * rho * U * U * L, "alpha": params.get("alpha", 0.0)}


if __name__ == "__main__":
    cond = sys.argv[1] if len(sys.argv) > 1 else "conditions.xml"
    r = reference(parameters(cond))
    print("  M = %.3f   Re = %.0f   alpha = %g deg" % (r["M"], r["Re"], r["alpha"]))
    print("  T_inf = %.2f K   c_inf = %.2f m/s   U_inf = %.3f m/s" % (r["T"], r["c"], r["U"]))
    print("  mu = %.6g Pa.s   nu = %.6g m2/s   q_inf = %.1f Pa" % (r["mu"], r["nu"], r["q"]))
