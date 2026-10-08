"""Overlay several burgers CSVs against the exact solution: python plot_compare.py out.png tag1 tag2 ..."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, 'scripts')
from plot_burgers import burgers_exact_sine, load
out, tags = sys.argv[1], sys.argv[2:]
fig, ax = plt.subplots(len(tags), 1, figsize=(7, 2.6*len(tags)), sharex=True, squeeze=False)
xe = np.linspace(-1, 1, 801)
for a, t in zip(ax[:,0], tags):
    x, u = load(f"results/{t}.csv")
    a.plot(xe, burgers_exact_sine(xe, 0.5), 'k--', lw=0.8); a.plot(x, u, lw=0.9); a.set_ylim(-1.6, 1.6); a.grid(alpha=.3)
    a.set_title(t, loc='left', fontsize=9)
fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)
