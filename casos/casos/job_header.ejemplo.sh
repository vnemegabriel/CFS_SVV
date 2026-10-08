#!/bin/bash
# Encabezado de cada trabajo en cola. `casos.py submit` reemplaza {id}, {procs} y {run}
# y agrega al final la línea que corre la corrida. Ejemplo para SLURM.
#SBATCH --job-name={id}
#SBATCH --ntasks={procs}
#SBATCH --time=72:00:00
#SBATCH --output={run}/job.out
