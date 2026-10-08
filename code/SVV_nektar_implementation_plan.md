# SVV en Nektar++ CompressibleFlowSolver — mapa de archivos

Qué hay que leer, entender y tocar. Sin justificaciones.

Raíz Nektar++ 5.9.0: `/opt/nektar++` (WSL). Rutas relativas a `library/` o `solvers/`.

---

## 1. Cadena de ejecución (leer en este orden)

Es el camino que recorre un paso de tiempo. Hay que poder seguirlo entero antes de tocar nada.

| #   | Archivo                                                                                  | Qué hace                                                                                                        |
| --- | ---------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| 1   | `solvers/CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.cpp` → `DoOdeRhs` | ensambla el RHS: advección + difusión + forcing                                                                 |
| 2   | mismo archivo → `DoOdeProjection`                                                        | proyecta la solución; acá se aplica el filtro exponencial si está activo                                        |
| 3   | `EquationSystems/NavierStokesCFE.cpp` → `v_DoDiffusion`                                  | término viscoso (LDGNS o InteriorPenalty) y llamada a viscosidad artificial                                     |
| 4   | `ArtificialDiffusion/ArtificialDiffusion.{h,cpp}`                                        | clase base + factory; dueña del operador LDG sobre variables conservadas                                        |
| 5   | `ArtificialDiffusion/NonSmoothShockCapture.cpp`                                          | ejemplo completo de plugin: define `mu(x)` con sensor de Persson                                                |
| 6   | `library/SolverUtils/Diffusion/DiffusionLDG.cpp`                                         | de dónde sale `qfield = ∇u`                                                                                     |
| 7   | `library/StdRegions/StdQuadExp.cpp` → `v_SVVLaplacianFilter`                             | el kernel SVV ya implementado (copiar este patrón)                                                              |
| 8   | `library/StdRegions/StdRegions.hpp`                                                      | enums `eFactorSVVCutoffRatio`, `eFactorSVVDiffCoeff`, `eFactorSVVPowerKerDiffCoeff`, `eFactorSVVDGKerDiffCoeff` |

Variantes por forma de elemento del punto 7 (necesarias para tri/hex): `StdTriExp.cpp`, `StdHexExp.cpp`, `StdTetExp.cpp`, `StdPrismExp.cpp`.

---

## 2. Conceptos de la biblioteca que hay que saber

Mínimo indispensable, en orden de dependencia:

1. **`ExpList` / `LocalRegions::Expansion` / `StdRegions::StdExpansion`** — jerarquía: lista global de expansiones → elemento físico → elemento de referencia. `m_fields[i]->GetExp(n)` da el elemento `n`.
2. **`FwdTrans` / `BwdTrans` / `PhysDeriv`** — pasaje coeficientes ↔ puntos físicos. Es la operación central del filtro modal.
3. **Bases ortogonales `eOrtho_A/B/C`** — la base modal donde el kernel SVV es diagonal. En triángulos/tetraedros hay coordenadas colapsadas; `v_SVVLaplacianFilter` ya resuelve eso.
4. **`StdMatrixKey` + `ConstFactorMap`** — cómo se pasan los parámetros SVV a los operadores existentes.
5. **`Array<OneD, ...>`** — contenedor de Nektar++; semántica de referencia, no de copia.
6. **Patrón factory (`NekFactory`)** — cómo se registra un plugin nuevo y cómo se lo selecciona desde el XML de sesión.
7. **Sesión XML: `SOLVERINFO` y `PARAMETERS`** — dónde se leen las claves nuevas (`m_session->DefinesSolverInfo(...)`, `m_session->LoadParameter(...)`).

---

## 3. Qué hay que tocar

| Archivo | Acción |
|---|---|
| `CompressibleFlowSolver/ArtificialDiffusion/SVVDiffusion.{h,cpp}` | **nuevo** — registrar en `ArtificialDiffusionFactory` como `"SVV"`; sobrescribir `GetFluxVector` con el kernel modal sobre `qfield`; leer `SVVDiffCoeff`, `SVVCutoffRatio`, tipo de kernel |
| `CompressibleFlowSolver/CMakeLists.txt` | agregar los fuentes nuevos |
| `EquationSystems/CompressibleFlowSystem.cpp` | instanciar el plugin desde su propia clave de sesión |
| `EquationSystems/EulerCFE.h`, `NavierStokesCFE.h` | `v_SupportsShockCaptType` si se usa la ruta de alias por `ShockCaptureType` |
| `CompressibleFlowSolver/Tests/` | tests de regresión propios: estado constante, transparencia Couette, Poiseuille MMS (ver `CASOS.md`, nivel de verificación) |

Lo que cambia respecto de `NonSmoothShockCapture`: en vez de `viscousTensor[j][i] = mu(x) * qfield[j][i]`, se hace por elemento `FwdTrans` ortogonal de `qfield` → escalar modo `k` por `ε·Q̂ₖ` → `BwdTrans`.

También hay que cubrir la ruta `DiffuseCoeffs` (espejo de `v_DoArtificialDiffusionCoeff`) para no romper implícito/ALE.

---

## 4. Claves de sesión que ya existen (baselines, sin tocar código)

- `SPECTRALHPDEALIASING` — sobreintegración
- `ExponentialFiltering` + `FilterAlpha` / `FilterExponent` / `FilterCutoff`
- `ShockCaptureType = NonSmooth` — AV laplaciana con sensor de Persson
- `ShockCaptureType = Physical` — AV física (requiere difusión IP)

---

## 5. Repos locales de apoyo

| Repo | Archivos relevantes |
|---|---|
| `TESIS/code/ndg_methods` | `Codes1D/Filter1D.m` (filtro modal), `CFD1D/Euler*` (Sod, Shu–Osher), `CFD2D/EulerShock2D.m` + `ForwardStepBC2D/IC2D.m`, `CurvedCNS2D.m` |
| `TESIS/code/ESDG` | `examples/IDP/dg1D_euler_shuosher.jl`, casos 2D CNS |
| `TESIS/code/3d_mortensen_dns` | `mortensens.py`, `tgv_final.h5` — referencia DNS para TGV |
