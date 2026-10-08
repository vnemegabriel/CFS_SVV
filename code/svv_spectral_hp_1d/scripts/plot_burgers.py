"""Reproduce K&S Fig. 6.27 from the CSV output of burgers_svv.

usage: python scripts/plot_burgers.py [results/burgers_nosvv.csv results/burgers_svv.csv]
"""
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def burgers_exact_sine(x, t):
    """u_t + u u_x = 0, u0 = -sin(pi x) on [-1,1]. Characteristics x = x0 + u0(x0) t.
    For t > 1/pi a standing shock sits at x = 0; each side is single-valued if
    we search x0 on the correct side of the origin."""
    u = np.zeros_like(x)
    for i, xi in enumerate(x):
        if abs(xi) < 1e-14:
            u[i] = 0.0
            continue
        a, b = (-1.0, 0.0) if xi < 0 else (0.0, 1.0)
        # x(x0) is monotone on the branch that reaches the shock first; bisection
        for _ in range(80):
            m = 0.5 * (a + b)
            if m - np.sin(np.pi * m) * t < xi:
                a = m
            else:
                b = m
        x0 = 0.5 * (a + b)
        u[i] = -np.sin(np.pi * x0)
    return u


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["x"], d["u"]


def load_spectrum(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    nel = int(d["elem"].max()) + 1
    nk = int(d["k"].max()) + 1
    return np.abs(d["coef"]).reshape(nel, nk)


def main():
    f_nosvv = sys.argv[1] if len(sys.argv) > 1 else "results/burgers_nosvv.csv"
    f_svv = sys.argv[2] if len(sys.argv) > 2 else "results/burgers_svv.csv"

    fig, ax = plt.subplots(2, 1, figsize=(6.5, 6.5), sharex=True)
    xe = np.linspace(-1, 1, 801)
    ue = burgers_exact_sine(xe, 0.5)
    for a, f, title in zip(ax, [f_nosvv, f_svv], ["(a) without SVV", "(b) with SVV"]):
        x, u = load(f)
        a.plot(xe, ue, "k--", lw=0.8, label="exact")
        a.plot(x, u, "-", lw=1.0, label="CG spectral/hp")
        a.set_ylim(-1.5, 1.5)
        a.set_ylabel("u")
        a.set_title(title, loc="left", fontsize=10)
        a.grid(alpha=0.3)
        for xv in np.linspace(-1, 1, 6):
            a.axvline(xv, color="gray", lw=0.4, ls=":")
    ax[0].legend(loc="upper right", fontsize=8)
    ax[1].set_xlabel("x")
    fig.suptitle("Inviscid Burgers, T = 0.5, 5 elements x 16 modes  (K&S Fig. 6.27)")
    fig.tight_layout()
    fig.savefig("results/fig_6_27_burgers.png", dpi=150)
    print("wrote results/fig_6_27_burgers.png")

    # modal spectra: what the kernel actually does
    try:
        s0 = load_spectrum(f_nosvv.replace(".csv", "_spectrum.csv"))
        s1 = load_spectrum(f_svv.replace(".csv", "_spectrum.csv"))
    except OSError:
        return
    nel, nk = s0.shape
    fig, ax = plt.subplots(1, nel, figsize=(3.0 * nel, 3.2), sharey=True)
    k = np.arange(nk)
    for e in range(nel):
        ax[e].semilogy(k, s0[e] + 1e-16, "o-", ms=3, label="no SVV")
        ax[e].semilogy(k, s1[e] + 1e-16, "s-", ms=3, label="SVV")
        ax[e].axvline(8, color="gray", ls=":", lw=0.8)
        ax[e].set_title(f"element {e}", fontsize=9)
        ax[e].set_xlabel("Legendre mode k")
        ax[e].grid(alpha=0.3)
    ax[0].set_ylabel("|u~_k|  (orthonormal coefficient)")
    ax[0].legend(fontsize=8)
    fig.suptitle("Per-element Legendre spectra at T=0.5 (dotted line: M_SVV = 8)")
    fig.tight_layout()
    fig.savefig("results/burgers_spectra.png", dpi=150)
    print("wrote results/burgers_spectra.png")


if __name__ == "__main__":
    main()
