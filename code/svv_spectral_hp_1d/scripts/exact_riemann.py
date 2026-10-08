"""Exact Riemann solver for the 1D Euler equations (Toro, ch. 4).

exact_riemann(rhoL,uL,pL, rhoR,uR,pR, gamma, x, t, x0) -> rho, u, p arrays.
"""
import numpy as np


def _f_k(p, rhok, pk, gamma):
    ck = np.sqrt(gamma * pk / rhok)
    if p > pk:  # shock
        A = 2.0 / ((gamma + 1) * rhok)
        B = (gamma - 1) / (gamma + 1) * pk
        f = (p - pk) * np.sqrt(A / (p + B))
        df = np.sqrt(A / (p + B)) * (1 - (p - pk) / (2 * (p + B)))
    else:  # rarefaction
        f = 2 * ck / (gamma - 1) * ((p / pk) ** ((gamma - 1) / (2 * gamma)) - 1)
        df = 1.0 / (rhok * ck) * (p / pk) ** (-(gamma + 1) / (2 * gamma))
    return f, df


def star_state(rhoL, uL, pL, rhoR, uR, pR, gamma):
    p = 0.5 * (pL + pR)
    for _ in range(100):
        fL, dfL = _f_k(p, rhoL, pL, gamma)
        fR, dfR = _f_k(p, rhoR, pR, gamma)
        g = fL + fR + (uR - uL)
        dp = -g / (dfL + dfR)
        p_new = max(p + dp, 1e-10)
        if abs(p_new - p) < 1e-12 * p:
            p = p_new
            break
        p = p_new
    fL, _ = _f_k(p, rhoL, pL, gamma)
    fR, _ = _f_k(p, rhoR, pR, gamma)
    u = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    return p, u


def exact_riemann(rhoL, uL, pL, rhoR, uR, pR, gamma, x, t, x0=0.0):
    x = np.asarray(x, dtype=float)
    ps, us = star_state(rhoL, uL, pL, rhoR, uR, pR, gamma)
    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)
    rho = np.empty_like(x)
    u = np.empty_like(x)
    p = np.empty_like(x)
    s = (x - x0) / max(t, 1e-300)
    gm1, gp1 = gamma - 1, gamma + 1
    for i, si in enumerate(s):
        if si < us:  # left of contact
            if ps > pL:  # left shock
                rhoSL = rhoL * ((ps / pL) + gm1 / gp1) / (gm1 / gp1 * (ps / pL) + 1)
                SL = uL - cL * np.sqrt(gp1 / (2 * gamma) * ps / pL + gm1 / (2 * gamma))
                if si < SL:
                    rho[i], u[i], p[i] = rhoL, uL, pL
                else:
                    rho[i], u[i], p[i] = rhoSL, us, ps
            else:  # left rarefaction
                rhoSL = rhoL * (ps / pL) ** (1 / gamma)
                cSL = cL * (ps / pL) ** (gm1 / (2 * gamma))
                SHL, STL = uL - cL, us - cSL
                if si < SHL:
                    rho[i], u[i], p[i] = rhoL, uL, pL
                elif si > STL:
                    rho[i], u[i], p[i] = rhoSL, us, ps
                else:
                    fac = 2 / gp1 + gm1 / (gp1 * cL) * (uL - si)
                    rho[i] = rhoL * fac ** (2 / gm1)
                    u[i] = 2 / gp1 * (cL + gm1 / 2 * uL + si)
                    p[i] = pL * fac ** (2 * gamma / gm1)
        else:  # right of contact
            if ps > pR:  # right shock
                rhoSR = rhoR * ((ps / pR) + gm1 / gp1) / (gm1 / gp1 * (ps / pR) + 1)
                SR = uR + cR * np.sqrt(gp1 / (2 * gamma) * ps / pR + gm1 / (2 * gamma))
                if si > SR:
                    rho[i], u[i], p[i] = rhoR, uR, pR
                else:
                    rho[i], u[i], p[i] = rhoSR, us, ps
            else:  # right rarefaction
                rhoSR = rhoR * (ps / pR) ** (1 / gamma)
                cSR = cR * (ps / pR) ** (gm1 / (2 * gamma))
                SHR, STR = uR + cR, us + cSR
                if si > SHR:
                    rho[i], u[i], p[i] = rhoR, uR, pR
                elif si < STR:
                    rho[i], u[i], p[i] = rhoSR, us, ps
                else:
                    fac = 2 / gp1 - gm1 / (gp1 * cR) * (uR - si)
                    rho[i] = rhoR * fac ** (2 / gm1)
                    u[i] = 2 / gp1 * (-cR + gm1 / 2 * uR + si)
                    p[i] = pR * fac ** (2 * gamma / gm1)
    return rho, u, p


def sod(x, t):
    return exact_riemann(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4, x, t, 0.5)


if __name__ == "__main__":
    xx = np.linspace(0, 1, 11)
    print(np.c_[xx, np.array(sod(xx, 0.2)).T])
