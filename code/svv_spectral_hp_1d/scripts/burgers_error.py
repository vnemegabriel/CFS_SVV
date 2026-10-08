"""L2 error of a burgers_svv CSV against the exact (pre- or post-shock) solution."""
import sys, numpy as np
sys.path.insert(0, 'scripts')
from plot_burgers import burgers_exact_sine, load
t = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
x, u = load(sys.argv[1])
ue = burgers_exact_sine(x, t)
print(f"{sys.argv[1]}: max|u-ue|={np.max(np.abs(u-ue)):.3e}  rms={np.sqrt(np.mean((u-ue)**2)):.3e}  max|u|={np.max(np.abs(u)):.3f}")
