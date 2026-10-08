"""Solver falso para las pruebas: imita la salida de Nektar++ y falla a pedido."""
import math
import sys
import time

mode = sys.argv[1]
steps, dt = (400 if mode == "slow" else 200), 1e-3
print("EquationType: Fake", flush=True)
forces = open("forces.fce", "w")
forces.write("# Forces\n#   Time  F1-press  F1-visc  F1-total  F2-press  F2-visc  F2-total\n")

for n in range(steps + 1):
    t = n * dt
    cd = 0.02 + 0.01 * (math.sin(60 * t) if mode == "noconv" else math.exp(-t / 0.02))
    if n % 5 == 0:
        forces.write("%10g 0 0 %.10g 0 0 0.3\n" % (t, cd))
        forces.flush()
    if n and n % 10 == 0:
        print("Steps: %d        Time: %g        CPU Time: 0.01s" % (n, t), flush=True)
    if n == 50:
        if mode == "nan":
            print("Fatal   : Level 0 assertion violation\nNaN found during time integration.", flush=True)
            sys.exit(1)
        if mode == "crash":
            sys.exit(3)
        if mode == "fatal0":
            print("Fatal   : Level 0 assertion violation\nsomething else", flush=True)
            sys.exit(0)
        if mode == "hang":
            time.sleep(3600)
    time.sleep(0.02 if mode == "slow" else 0.001)
print("Time-integration  : 2s", flush=True)
