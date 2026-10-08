"""Plot rho, u, p from cns1d_svv runs, with exact solution overlays when known.

usage: python scripts/plot_cns.py --case sod results/sod_nosvv.csv results/sod_svv.csv [--out name]
       python scripts/plot_cns.py --case entropywave results/entropywave_svv.csv
"""
import argparse
import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from exact_riemann import sod  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--case", default="sod")
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    T = a.T if a.T is not None else {"sod": 0.2, "entropywave": 2.0, "shuosher": 1.8}[a.case]
    fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)

    exact = None
    if a.case == "sod":
        xe = np.linspace(0, 1, 1001)
        exact = (xe,) + tuple(sod(xe, T))
    elif a.case == "entropywave":
        xe = np.linspace(-1, 1, 801)
        exact = (xe, 1 + 0.2 * np.sin(np.pi * (xe - T)), np.ones_like(xe), np.ones_like(xe))

    if exact is not None:
        for k, name in enumerate(["rho", "u", "p"]):
            ax[k].plot(exact[0], exact[k + 1], "k--", lw=0.9, label="exact")

    for f in a.files:
        d = np.genfromtxt(f, delimiter=",", names=True)
        label = os.path.basename(f).replace(".csv", "")
        for k, name in enumerate(["rho", "u", "p"]):
            ax[k].plot(d["x"], d[name], "-", lw=1.0, label=label)

    for k, name in enumerate(["rho", "u", "p"]):
        ax[k].set_ylabel(name)
        ax[k].grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[2].set_xlabel("x")
    fig.suptitle(f"case = {a.case}, T = {T}")
    fig.tight_layout()
    out = a.out or f"results/cns_{a.case}.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
