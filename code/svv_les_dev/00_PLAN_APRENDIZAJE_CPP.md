# Plan de aprendizaje de C++ para implementar SVV y turbulencia en Nektar++

**Punto de partida asumido:** C++ casi nulo, con background en Python/MATLAB.
**Objetivo:** poder escribir, compilar y testear clases nuevas dentro de
`solvers/CompressibleFlowSolver/`.

---

## 0. Diagnóstico del repo (hecho, verificado por grep)

Antes de estudiar nada conviene saber exactamente qué existe y qué no.
Esto cambia radicalmente el alcance del trabajo.

| Qué | Dónde está en tu copia de Nektar++ | Estado |
|---|---|---|
| Kernel SVV bajo nivel (por elemento) | `library/StdRegions/StdExpansion.{h,cpp}` → `v_SVVLaplacianFilter`, `v_ExponentialFilter`, y su implementación en `StdQuadExp`, `StdHexExp`, `StdTriExp`, `StdTetExp`, `StdPrismExp`, `StdPyrExp`, `StdSegExp` | **Ya existe** |
| Bucle sobre elementos | `library/MultiRegions/ExpList.cpp:2352` → `ExpList::ExponentialFilter` | **Ya existe** |
| Orquestación SVV (parámetros, kernels exp/power/DG, coeficiente variable) | `solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp` → `SetUpSVV()`, `SVVVarDiffCoeff()`, `AppendSVVFactors()` (~líneas 993–1245) | **Ya existe, pero solo para el solver INCOMPRESIBLE** |
| SVV en `CompressibleFlowSolver` | — | **No existe. Cero ocurrencias.** |
| Modelo de turbulencia (LES/RANS) en cualquier solver | — | **No existe.** Solo hay captura de choque (`ArtificialDiffusion/`) y generación de turbulencia sintética de entrada (`Forcing/ForcingCFSSyntheticEddy.cpp`) |

**Conclusión operativa:** tus dos tareas no son "configurar" nada. Son
**escribir clases C++ nuevas** que se enchufan en un framework existente.
Pero — y esto es la buena noticia — el framework te obliga a un molde muy
estrecho y repetitivo. No necesitás C++ "en general". Necesitás **un
subconjunto pequeño y muy identificable**.

---

## 1. El subconjunto de C++ que realmente necesitás

Abrí `solvers/CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.h`.
Son 135 líneas y contienen **casi todo** el C++ que vas a escribir en la tesis.
Esa es literalmente tu piedra de Rosetta.

Lista cerrada, en orden de aparición en ese archivo:

1. **Modelo de compilación**: `#ifndef`/`#define`/`#endif` (header guards),
   `#include`, la diferencia entre `.h` (declaración) y `.cpp` (definición),
   y por qué existe el linker.
2. **`namespace`** — todo Nektar++ vive en `namespace Nektar`.
3. **`typedef` / `using`** para dar nombre corto a tipos largos.
4. **`std::shared_ptr<T>`** — puntero con conteo de referencias. Nektar++ los
   usa en todos lados y casi nunca vas a ver `new`/`delete`.
5. **Clases**: `public` / `protected` / `private`, miembros de datos (prefijo
   `m_` por convención), constructor con **lista de inicialización**
   (`: m_session(pSession), m_fields(pFields)`).
6. **Referencias `&` y `const`**: la firma
   `const Array<OneD, const Array<OneD, NekDouble>> &inarray` no es magia — es
   "referencia constante a un array constante de arrays constantes de doubles".
   Entender *por qué* se escribe así (evitar copias + garantizar no-modificación)
   es el 20% del esfuerzo que da el 80% de la comprensión.
7. **Herencia y polimorfismo**: `class NonSmoothShockCapture : public ArtificialDiffusion`,
   `virtual`, `override`, `= 0` (función virtual pura → clase abstracta),
   destructor virtual.
8. **El idioma NVI (Non-Virtual Interface)**: método público `DoArtificialDiffusion()`
   que solo llama al protegido `v_DoArtificialDiffusion()`. **Todo Nektar++ está
   escrito así.** La `v_` es la que redefinís vos.
9. **Templates como consumidora, no como autora**: vas a *usar*
   `Array<OneD, NekDouble>` y `MemoryManager<X>::AllocateSharedPtr(...)`.
   **No** necesitás saber escribir templates.
10. **Miembros `static`**: el `static std::string className;` y el
    `static ... create(...)` que hacen funcionar el patrón factory.
11. **Punteros a función miembro**: aparece una sola vez, en
    `m_diffusion->SetFluxVector(&ArtificialDiffusion::GetFluxVector, this);`.
    Con reconocerlo alcanza.

### Lo que NO necesitás estudiar (y te va a tentar)

- Metaprogramación con templates, SFINAE, `constexpr` avanzado
- Semántica de movimiento (`std::move`, rvalue references) en profundidad
- Diseño de excepciones (Nektar++ usa sus propias macros `ASSERTL0/1`)
- Concurrencia / `std::thread` (el paralelismo acá es MPI, y ya está resuelto)
- La STL en profundidad (`<algorithm>`, iteradores exóticos). Nektar++ usa
  `Array<OneD,>` y `Vmath::` para los bucles numéricos, no `std::vector` +
  `std::transform`.

Esto reduce el temario de "6 meses" a **3–4 semanas de estudio dirigido**.

---

## 2. Fases

### Fase 0 — C++ mínimo viable (≈2 semanas, 1–2 h/día)

**Recurso principal:** [learncpp.com](https://www.learncpp.com/) — gratuito,
secuencial, con ejercicios. Es el mejor recurso para tu perfil.

Capítulos a hacer, **y parar ahí**:

- Cap. 1–2: fundamentos, funciones, archivos múltiples, header guards
- Cap. 4–5: tipos fundamentales, `const`
- Cap. 6: ámbito, `namespace`, duración de almacenamiento, `static`
- Cap. 9: referencias `&` y punteros `*`, `const` con referencias ← **crítico**
- Cap. 12–13: funciones con referencias, valores de retorno
- Cap. 14–15: clases, constructores, listas de inicialización, `this`
- Cap. 24: herencia
- Cap. 25: funciones virtuales, `override`, clases abstractas, destructores virtuales ← **crítico**
- Cap. 22: smart pointers (`std::shared_ptr`)

Salteá: sobrecarga de operadores en detalle, templates (cap. 11 y 26 — leelos
rápido solo para reconocer la sintaxis), lambdas, excepciones, contenedores STL.

**Herramienta de apoyo:** [Compiler Explorer](https://godbolt.org/) para probar
fragmentos de 10 líneas sin montar un proyecto.

**Criterio de aprobación de la fase:** podés explicar en voz alta, sin mirar,
qué hace cada línea de `ArtificialDiffusion.h`.

---

### Fase 1 — El ciclo compilar/editar/probar (≈1 semana)

Esta fase es **más importante que la Fase 0** y la gente la saltea. Si no tenés
un ciclo rápido de "cambio algo → compilo → veo el efecto", no vas a poder
aprender el codebase, porque el codebase se aprende rompiéndolo.

1. Compilar Nektar++ desde fuente en tu WSL, con `CMAKE_BUILD_TYPE=Debug` en
   un build separado (dejá también un build `Release` para correr casos).
   Activá `NEKTAR_SOLVER_COMPRESSIBLE_FLOW=ON` y `NEKTAR_BUILD_UNIT_TESTS=ON`.
2. Correr un caso existente de punta a punta:
   `solvers/CompressibleFlowSolver/Tests/` tiene decenas. Empezá por uno
   subsónico chico (`CylinderSubsonic_*`).
3. **Ejercicio deliberado:** meté un `std::cout` dentro de
   `NonSmoothShockCapture::v_GetArtificialViscosity`, recompilá, corré un caso
   con captura de choque, y verificá que se imprime. Cronometrá cuánto tarda el
   ciclo completo. Si tarda más de ~5 minutos, aprendé a compilar solo el
   target del solver (`make CompressibleFlowSolver -j`) en vez de todo.
4. Aprendé `ccache` — te va a ahorrar semanas acumuladas.
5. Aprendé lo mínimo de `gdb`: `break`, `run`, `bt`, `print`. O configurá VS Code
   con la extensión de C++ y el `launch.json` apuntando al ejecutable del solver.

**Criterio de aprobación:** podés poner un breakpoint en una función del solver
y ver el stack de llamadas que llevó hasta ahí.

---

### Fase 2 — El dialecto Nektar++ (≈1–2 semanas)

Acá no estudiás C++ estándar, estudiás las convenciones de la casa. Leé en
**este orden exacto**:

| # | Archivo | Qué extraer |
|---|---|---|
| 1 | `library/LibUtilities/BasicUtils/SharedArray.hpp` | Qué es `Array<OneD, T>`: array con conteo de referencias. Cómo se construye (`Array<OneD,NekDouble>(n, 0.0)`), qué significa `Array<OneD, Array<OneD, NekDouble>>` (array de campos), y por qué copiar uno **no** copia los datos |
| 2 | `library/LibUtilities/BasicUtils/NekFactory.hpp` | El patrón factory: cómo un string del XML de sesión se convierte en un objeto. Es el mecanismo por el que vas a "enchufar" tus clases |
| 3 | `library/LibUtilities/BasicUtils/VmathArray.hpp` y `Vmath.hpp` | `Vmath::Vadd`, `Vmath::Vmul`, `Vmath::Smul`, `Vmath::Vvtvp`. **Toda la aritmética vectorizada se escribe con esto, no con bucles a mano** |
| 4 | `solvers/CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.{h,cpp}` | La piedra de Rosetta: factory + NVI + constructor protegido |
| 5 | `solvers/CompressibleFlowSolver/ArtificialDiffusion/NonSmoothShockCapture.{h,cpp}` | Cómo se ve una clase derivada concreta: `create()`, `className`, `v_...() override` |
| 6 | `solvers/CompressibleFlowSolver/CMakeLists.txt` | Cómo se agregan archivos nuevos al build (líneas 11–12 y 71–72 son el par `.cpp`/`.h` de ArtificialDiffusion) |
| 7 | `solvers/CompressibleFlowSolver/EquationSystems/NavierStokesCFE.{h,cpp}` | Dónde vive la física viscosa. Buscá `GetViscousFluxVector`, `GetViscosityAndThermalCondFromTemp`, `v_DoDiffusion`. **Acá es donde se engancha la turbulencia** |
| 8 | `solvers/CompressibleFlowSolver/Misc/VariableConverter.h` | Cómo se pasa de variables conservativas (ρ, ρu, ρE) a temperatura, presión, velocidad. Lo vas a necesitar en ambas tareas |
| 9 | `solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp` (`SetUpSVV`, `SVVVarDiffCoeff`, `AppendSVVFactors`) | **La referencia canónica de SVV en este código.** No la copies literal (el solver incompresible aplica SVV dentro del operador de Helmholtz; el compresible es explícito y necesita otro enganche), pero copiá los nombres de parámetros y la lógica de los kernels |

**Ejercicio de la fase:** escribí a mano, en papel, el diagrama de llamadas
desde `CompressibleFlowSystem::DoOdeRhs` hasta el bucle sobre elementos. Vas a
usarlo todo el resto de la tesis.

---

### Fase 3 — Escribir el código (el resto)

En este punto ya no es aprendizaje de C++, es trabajo de tesis. El esqueleto
está en `skeleton/` junto a este archivo, con comentarios línea por línea.

Orden sugerido:

1. **SVV primero.** Es más acotado, el kernel bajo nivel ya existe en
   `StdRegions`, y no toca la física — es un operador numérico. Te da la victoria
   temprana que necesitás para no frustrarte.
2. **Turbulencia después.** Requiere entender el tensor de deformación, el
   filtrado en flujo compresible (Favre), y decisiones físicas (¿Smagorinsky?
   ¿Vreman? ¿WALE? ¿corrección de compresibilidad de Yoshizawa?). El C++ es más
   fácil que en SVV; la física es mucho más difícil.

---

## 3. Nota sobre la parte sub → supersónica

El pedido "modelen turbulencia adecuadamente para casos sub a supersónicos" es
donde está el riesgo real de la tesis, y **no es un problema de C++**. Tres cosas
a tener presentes desde el diseño de las clases:

- **Interacción SVV ↔ captura de choque.** En supersónico ya hay
  `ArtificialDiffusion` metiendo viscosidad artificial cerca de los choques. Si
  además metés SVV y un modelo LES, tenés tres fuentes de disipación numérica
  compitiendo. El diseño debe permitir **activarlas y medirlas por separado**
  (por eso en el esqueleto cada una expone su propio campo de viscosidad).
- **Sensor de choque vs. sensor de turbulencia.** Un Smagorinsky estándar no
  distingue un choque de una capa de corte turbulenta, y le mete viscosidad
  turbulenta enorme al choque. Vreman y WALE se comportan mejor; considerá
  también un limitador basado en dilatación (`∇·u < 0` fuerte ⇒ choque).
- **Promediado de Favre.** En compresible el filtrado es ponderado por densidad.
  Si escribís las fórmulas incompresibles con ρ constante, el modelo va a estar
  mal en supersónico aunque compile perfecto.

---

## 4. Estrategia de verificación (aplicable desde el día 1)

Tenés ya en `code/3d_mortensen_dns/` un caso Taylor–Green con
`energy_spectrum.png`. **Ese es tu banco de pruebas ideal para SVV y LES**,
porque el efecto de ambos es precisamente sobre el extremo de altas frecuencias
del espectro de energía.

Escalera de tests, de más barata a más cara:

| Nivel | Qué verifica | Cómo |
|---|---|---|
| 1. Unit test | Que el kernel SVV hace lo que decís | `library/UnitTests/` — proyectá un modo conocido, aplicá el filtro, comprobá la amplitud resultante contra la fórmula analítica del kernel |
| 2. Test de regresión | Que no rompiste nada | `.tst` en `solvers/CompressibleFlowSolver/Tests/` + `ctest`. Copiá un `.tst` existente y adaptalo |
| 3. Consistencia | Que con coeficiente SVV → 0 recuperás el solver original bit a bit | Correr el mismo caso con `SVVDiffCoeff = 0` y diferenciar los `.fld` |
| 4. Orden de convergencia | Que no destruiste la precisión espectral | Caso con solución analítica (vórtice isentrópico — ya está en `IsentropicVortexBC`), barrer P = 2..8, verificar pendiente |
| 5. Física, subsónico | Que el modelo hace lo correcto | Taylor–Green Re=1600: comparar disipación de energía cinética vs. DNS de referencia; espectro con pendiente −5/3 |
| 6. Física, supersónico | Que sobrevive a los choques | Turbulencia isotrópica compresible decayente (Samtaney et al.) barriendo Mach turbulento; y un caso con choque para verificar que el modelo no arruina el salto de Rankine–Hugoniot |

Los niveles 1–3 son los que te dan confianza para seguir programando. Los 4–6
son los que van al escrito.

---

## 5. Ritmo realista

| Semanas | Qué |
|---|---|
| 1–2 | Fase 0 (C++ mínimo) |
| 3 | Fase 1 (build + debug + ciclo rápido) |
| 4–5 | Fase 2 (dialecto Nektar++, lectura dirigida) |
| 6–8 | SVV: implementación + tests 1–4 |
| 9–14 | Turbulencia: implementación + calibración + tests 5–6 |

Si en la semana 3 todavía no compilaste el solver, parate y resolvé eso antes de
seguir con teoría: es el cuello de botella real.
