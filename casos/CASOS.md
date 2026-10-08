# Casos
## 1. Verificación — unit tests de Nektar++

`ctest` sobre `CompressibleFlowSolver` con SVV compilado. Si pasan todos, el operador no rompió nada.
Criterio binario, sin interpretación.

```bash
cd ~/nektar-master/build
ctest -R Compressible -j8
```

## 2. Calibración — Shu–Osher y tubo de choque

Acá se fijan `μ_SVV`, `r_cut`, núcleo, variable del sensor y umbral. Con barridos.

| Caso | Sesiones | Qué aporta |
|---|---|---|
| Shu–Osher 1D | 50 elem P5, 100 elem P5 y P8 (Mateo-Gabín et al. 2022 §6.2) | choque fuerte + estructura fina; referencia ENO-RF-S-3, 1600 gdl |
| Tubo de choque 2D | `ShockTube_2D_mixedMesh_*`, 14 sesiones | triángulos; 8 solvers de Riemann sobre configuración idéntica separan disipación de upwinding de la del filtro |

Los tres núcleos están en `StdRegions`; la orquestación hay que portarla del solver incompresible.

## 3. Costo–precisión: NACA 0012

El resultado central.

|                           | Mallas | Órdenes | Mach                |
| ------------------------- | ------ | ------- | ------------------- |
| Nektar++                  | 3      | 3       | 0.3, 0.8, 1.3, 2, 3 |
| OpenFOAM `rhoCentralFoam` | 3      | —       | ídem                |

45 corridas por estabilización (DG puro, viscosidad artificial, SVV) y 15 de OpenFOAM.

- `Re` crece con el Mach: `ρ`, `p` y `μ` fijos, `cInf` constante, el Mach escala la velocidad.
- Eje común: grados de libertad = elementos·(P+1)² en Nektar++, celdas en OpenFOAM.
- Un solo ángulo de ataque. Un solo radio de campo lejano.
- Cantidad de interés: `Cd` primaria, `Cl` secundaria.
- Criterio de parada idéntico en ambos códigos, sobre `Cd`.

## 4. Aplicación: Aconcagua

Contra OpenFOAM y contra telemetría. `Ma 0 a 1.8`, misma malla en todos los puntos.
Se mide `C_D`, distribución de presión, posición del choque de proa y presión de base.

