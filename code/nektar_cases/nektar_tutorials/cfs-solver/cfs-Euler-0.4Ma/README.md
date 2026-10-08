# Nektar++ sequencing

geometry file: `naca.xml`
condition file: `session_naca.xml`

## run sim
```bash
mpirun -np 4 CompressibleFlowSolver naca.xml sessiom_naca.xml
```

## postprocessing

```bash
python3 postNek.py --mesh naca.xml --session session_naca.xml #ver script para modificaciones
python3 viewNek.py vtu/naca.pvd --field u
```
