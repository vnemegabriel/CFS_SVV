	# Estructura de la tesis

> *¿De qué forma se puede lograr, con un solver basado en métodos espectrales, la misma robustez y precisión que un solver FVM estándar, a menor costo de cómputo?*
> 
> La propuesta de respuesta es la estructuración de un flujo de trabajo el operador de **viscosidad espectral desvaneciente (SVV)** como mecanismo de estabilización y robustez; la verificación es la curva costo–error contra OpenFOAM
> en el rango M = 0.3–3, con el cohete **Aconcagua** como caso de aplicación y su telemetría de
> vuelo como referencia experimental.

---

## 0. Orden global de capítulos

```latex
% ---------------- Parte 1: marco argumental ----------------
\include{INTRODUCCION}          % contexto → problema → aporte
\include{ESTADO_DEL_ARTE}       % termina identificando EL hueco
\include{OBJETIVOS_Y_ALCANCE}   % hueco → pregunta → qué NO se hace
% ---------------- Parte 2: teoría --------------------------
\include{MARCO_TEORICO}         % NSE compresibles + DGSEM spectral/hp
\include{TURBULENCIA_Y_ESCALAS} % qué se resuelve, qué se modela, qué se disipa
\include{FENOMENOS}             % qué rompe el esquema y por qué
\include{FILTRO_SVV}            % el operador: derivación, kernel, sensor, parámetros
% ---------------- Parte 3: implementación ------------------
\include{BIBLIOTECA_NEKTAR}     % arquitectura, inserción, verificación de código
\include{CASOS}                 % diseño de los experimentos numéricos
% ---------------- Parte 4: cierre --------------------------
\include{RESULTADOS}
\include{DISCUSION}
\include{CONCLUSION}
```

---

# PARTE 1 — MARCO ARGUMENTAL

## Capítulo 1 — `INTRODUCCION`

### 1.1 El CFD en el diseño aerodinámico y su límite actual

- El CFD pasó de herramienta de verificación a herramienta de diseño: las decisiones de peso se
  toman sobre simulaciones, no sobre túnel.
- El costo de fidelidad sigue rendimientos decrecientes. Fuente dura para cuantificarlo:
  Slotnick et al. (2014), *CFD Vision 2030*, NASA/CR-2014-218178 ** **, que estima el costo de
  LES con modelo de pared escalando como `Re^1.3` en FLOPs y de LES resuelta en pared como
  `Re^2.5`, y concluye que *"LES methods will remain impractical for various important
  applications for the foreseeable future"*.
- El mismo informe reporta que la LES con modelo de pared es viable en 24 h a `Re ≈ 10^6` sobre
  geometrías de relación de aspecto unitaria con ~6 PFLOP/s, y que para relaciones de aspecto
  realistas de flujo externo el costo es un orden de magnitud mayor, *"out of the reach of even
  2030 leadership HPC machines"* ****.
- Consecuencia distributiva: la simulación de alta fidelidad viable en producción queda reservada
  a grandes grupos e industria pesada. Enunciar el problema en términos de **acceso**: la pregunta
  no es "¿cuál es el método más exacto?" sino "¿cuánta precisión se compra por hora de núcleo
  disponible?".
- Aterrizar en el caso local: un equipo universitario que diseña un cohete supersónico, tiene
  telemetría de vuelo y necesita coeficientes aerodinámicos confiables sin acceso a cómputo
  masivo. Ese caso reaparece en la Parte 3 — plantarlo acá.

### 1.2 Métodos de alto orden: promesa y fricción

- FVM de segundo orden como estándar de facto: robustez, geometrías complejas, madurez de
  software, tolerancia a mallas malas.
- La promesa del alto orden: menos grados de libertad para una tolerancia dada, y por lo tanto
  menor costo. Evidencia citable: Wang et al. (2013), *IJNMF* 72:811–845 ****, primer taller
  internacional de métodos de alto orden, que reporta que *"High-order methods demonstrated better
  performance than the 2nd-order finite volume method for both steady and unsteady problems based
  on error versus cost"*.
- **Y su límite documentado, que hay que citar en la misma página** — del mismo trabajo ****:
  *"For problems with non-smooth solutions or geometries, high-order methods cannot achieve
  high-order accuracy, as expected. There is not enough evidence to state whether high-order
  methods perform better or worse than low-order methods. More comparisons are necessary."*
  Ese hueco —soluciones no suaves, es decir choques— es exactamente el rango de esta tesis.
- La fricción principal: **fragilidad**. Baja disipación significa que nada amortigua ni la
  energía que se acumula en los modos altos ni las oscilaciones de Gibbs. Evidencia directa y
  reciente: Rubio et al. (2025), arXiv:2512.04574 ****, reportan que en DG estándar sin forma
  split, sobre el vórtice de Taylor–Green a Re = 1600, *"the central configuration becomes
  unstable around t/t_c ≈ 4, while the Roe variant remains stable for a longer period but fails
  near t/t_c ≈ 7  The Riemann solver alone does not provide sufficient dissipation to stabilize
  the simulation"*. **Ésta es la configuración de esta tesis: DG estándar. La cita es la
  motivación más directa que existe.**
- Fricciones secundarias, mencionar sin desarrollar: generación de mallas curvadas —Wang et al.
  (2013) la listan como *pacing item* número uno ****—, restricción de paso temporal explícito,
  curva de aprendizaje del software.

### 1.3 La búsqueda de la precisión barata

- Qué se hace hoy para domar la fragilidad
	- filtrado modal
	- viscosidad artificial con sensor
	- limitadores
	- esquemas entrópicamente estables con forma split y suma por partes
	- captura sub-celda.
- El patrón que se repite: **un mecanismo distinto para cada patología**, con proliferación de
  constantes de ajuste. Evidencia en la fuente original del mecanismo que usa Nektar++ —
  Persson y Peraire (2006), AIAA 2006-112 ****: *"κ is chosen empirically"*; *"for high Mach
  numbers, the value of the viscosity coefficient needs to be increased to maintain stability and
  this results in wider shocks"*; *"the value of the viscosity coefficient has been tuned for each
  Mach number"*.
- Dónde funciona hoy el alto orden: aeroacústica, LES de baja velocidad, turbomaquinaria. El
  compresible externo con choques queda rezagado.
- La alternativa de esta tesis: un solo operador con teoría de convergencia.

### 1.4 Aporte de este trabajo

- Enunciado en una frase: se implementa un operador de viscosidad espectral desvaneciente
  modulado por sensor modal en el `CompressibleFlowSolver` de Nektar++, se caracteriza su comportamiento paramétrico en M = 0.3–3, y se mide el costo a precisión fija del esquema resultante frente a OpenFOAM sobre el cohete Aconcagua.
- Aportes, en orden de peso:
	1. Implementación del operador SVV en el solver compresible de Nektar++, cerrando una brecha verificada en el código: los núcleos SVV existen en `library/StdRegions` para todas las formas de elemento, pero son inalcanzables desde el solver compresible (ver 2.5 y 8.2).
	2. Caracterización sistemática de los parámetros del operador frente al número de Mach.
	3. Comparación contra los mecanismos de estabilización ya disponibles en el propio solver.
	4. Curva costo–error contra OpenFOAM sobre geometría real, con telemetría de vuelo como referencia.
- Declarar de entrada lo que **no** es aporte: la formulación del operador se adopta de la
  literatura; no se desarrolla teoría nueva, y no se demuestra estabilidad entrópica (ver 3.5 y
  Anexo A.1).

---

## Capítulo 2 — `ESTADO_DEL_ARTE`

### 2.1 Métodos espectrales y de alto orden para flujo compresible

- Linaje: métodos espectrales globales → limitación geométrica → elementos espectrales y
  spectral/hp (Karniadakis y Sherwin, 2005 ****) → DG nodal (Hesthaven y Warburton, 2008
  ****) y su variante colocada DGSEM.
- Base modal jerárquica frente a nodal; por qué Nektar++ adopta la modal jerárquica y qué implica
  para SVV: el filtro es diagonal en la base modal ortogonal, que es exactamente lo que el
  operador necesita.
- Estado de madurez del software: Nektar++ (Cantwell et al., 2015, *CPC* 192:205–219 ****;
  Moxey et al., 2020, *CPC* 249:107110 ****), Nek5000, FLEXI, HORSES3D, Trixi.jl.
- Observación estructural que después importa: los códigos DGSEM con forma split (HORSES3D,
  FLUXO, Trixi.jl) son tensoriales hexaédricos por diseño, mientras Nektar++ es deliberadamente
  agnóstico en base y admite tri/tet/prisma/pirámide. Esa diferencia de diseño es la razón por la
  que la ruta entrópicamente estable no está disponible acá (ver Anexo A.1).

### 2.2 Turbulencia compresible: por qué no es turbulencia incompresible con densidad variable

*(Este capítulo del estado del arte se limita a ubicar la literatura. El desarrollo va en el cap. 5.)*

- Formalismo de LES compresible: filtrado de Favre y proliferación de términos submalla en la
  ecuación de energía (Erlebacher et al., 1992, *JFM* 238:155–185 ****; Garnier, Adams y
  Sagaut, 2009 ****).
- Número de Mach turbulento y disipación dilatacional (Sarkar et al., 1991, *JFM* 227:473–493
  ****; Zeman, 1990 ****; Blaisdell, Mansour y Reynolds, 1993, *JFM* 256:443–485
  ****). Revisión: Lele (1994), *Annu. Rev. Fluid Mech.* 26:211–254 ****.
- **El dato que ordena todo el capítulo 5** — Duan, Beekman y Martín (2011), *JFM* 672:245–267
  ****: en capas límite turbulentas de M = 0.3 a 12, *"the explicit dilatation terms such as
  pressure dilatation and dilatational dissipation remain small for the present Mach number
  range"* y las relaciones de escalado de Morkovin siguen valiendo. Es decir: para un cuerpo
  esbelto con capa límite adherida en el rango de esta tesis, **la compresibilidad que importa es
  la del flujo medio, no la de la turbulencia**.

### 2.3 Viscosidad espectral desvaneciente: linaje teórico

- Tadmor (1989), *SIAM J. Numer. Anal.* 26:30–44 ****: convergencia de métodos espectrales para
  leyes de conservación no lineales mediante viscosidad espectral. Es el único fundamento de
  convergencia que tiene cualquiera de los mecanismos considerados. Alcance real: leyes de
  conservación escalares 1D en base de Fourier — decirlo.
- Maday, Ould Kaber y Tadmor (1993), *SIAM J. Numer. Anal.* 30:321–342 ****: núcleo de
  Legendre, el aplicable a discretizaciones por elementos. Alcance: escalar 1D y dinámica de gases
  1D.
- Kirby y Sherwin (2006), *CMAME* 195:3128–3144 ****: adaptación al marco spectral/hp, por los
  propios autores de Nektar++. **Incompresible.** Y notablemente reservados: su conclusión es que
  SVV *"provides another computational tool for stabilisation"*, no que sea un modelo LES. También
  advierten que *"SVV is not necessarily a means of fixing all errors associated with problems
  such as consistent-integration"*.
- Moura, Sherwin y Peiró (2016), *JCP* 307:401–422 ****: análisis de autosoluciones que provee
  criterios analíticos de parámetros e introduce el núcleo de ley de potencias. Resultado citable:
  las características de disipación de CG con SVV son *"broadly equivalent to those obtained
  through upwinding in the discontinuous Galerkin (DG) scheme"*. También corrige un defecto del
  escalado tradicional, que hace depender el operador de un número de Péclet local.
- **Toda esta línea es incompresible o escalar lineal.** Karamanos y Karniadakis (2000), Kirby y
  Karniadakis (2002), Pasquetti (2006), Kirby y Sherwin (2006): incompresible sin excepción
  ****. Decirlo explícitamente es honestidad barata y es además parte del argumento de novedad.

### 2.4 SVV en flujo compresible y en DG

- Manzanero, Ferrer, Rubio y Valero (2020), *Computers & Fluids* 200:104440 ****: primer
  modelo de turbulencia DG-SVV; amplitud dada por Smagorinsky (`μ_SVV = C_S² Δ² |S|`, `C_S = 0.2`)
  y filtrada por el núcleo de ley de potencias de Moura, con óptimo `P_SVV = 0.1`. Formalmente
  compresible. **Pero**: *"In all cases, the Mach number is kept to 0.1 such that compressible
  effects are negligible"* ****.
- Mateo-Gabín, Manzanero y Valero (2022), *JCP* 471:111618 ****: SVV entrópicamente estable en
  DG; casos Shu–Osher, escalón supersónico M = 3 y vórtice de Taylor–Green. El TGV que corren es
  el **inviscido**, a M ≈ 0.09. **Corrección de bibliografía: son tres autores —Mateo-Gabín,
  Manzanero y Valero—, no cuatro; el archivo local `Gabin2008_...pdf` tiene mal el año y la
  entrada previa de esta tesis tenía mal la lista de autores.**
- Dato que hay que citar: los propios autores reportan **** que *"SVV filtering is not enough
  to control the oscillations near the strong shock"* — de ahí la necesidad del sensor.
- El hueco que queda: **toda la literatura de SVV compresible corre a Mach efectivamente
  incompresible.** La formulación se validó sobre ecuaciones compresibles; la física, no.
- Alternativas contemporáneas al mismo problema, para el contraste: penalización de saltos de
  gradiente (Moura et al., 2022, *CMAME* 388:114200 ****; Kou, Mariño y Ferrer, 2022
  ****), LES modal multiescala en DG (Chapelier, de la Llave Plata y Lamballais, 2016, *CMAME*
  307:275–299 **[F]**).

### 2.5 Estabilización disponible hoy en Nektar++ y sus límites

*(Sección corta acá; el inventario completo va en el cap. 8. Su función es dejar planteado el hueco.)*

- Verificado sobre el árbol de fuentes de Nektar++ v5.10.0 ****: el `CompressibleFlowSolver`
  estabiliza por cuatro vías — disipación de upwind del flujo de Riemann, viscosidad artificial
  laplaciana con sensor de Persson–Peraire (`ShockCaptureType = NonSmooth`), viscosidad artificial
  física (`Physical`), filtro modal exponencial y desaliasing por sobre-integración.
- `grep -rn "SVV" solvers/CompressibleFlowSolver/` devuelve **cero coincidencias** ****. Los
  núcleos SVV sí existen en `library/StdRegions` para todas las formas de elemento, pero son
  alcanzables sólo a través del operador de Helmholtz del solver incompresible; el solver
  compresible nunca ensambla un Helmholtz. **SVV no está sin usar: está inalcanzable.**
- La asimetría es citable: Moxey et al. (2020) **** presentan SVV entre las respuestas de
  Nektar++ a la falta de robustez —*"efficient dealiasing techniques and spectral vanishing
  viscosity, that have proved invaluable for particularly challenging applications"*— y esa
  afirmación es cierta sólo del lado incompresible.
- Verificado también ****: Nektar++ no tiene ninguna discretización entrópicamente estable,
  ni forma split, ni suma por partes, ni flujos de dos puntos, ni limitador de positividad.

### 2.6 Comparaciones de costo entre alto orden y segundo orden

- La metodología aceptada. Wang et al. (2013) ****, verbatim: *"We believe that this
  [error versus cost] is a more objective way to measure different methods than to compare the
  error on a given mesh or even with the same degrees of freedom."*
- Unidad de costo estandarizada: *work units* normalizadas por TauBench, excluyendo inicialización,
  post-proceso y E/S (guías HiOCFD4 y HiOCFD5 ****). Regla de prudencia del mismo trabajo
  ****: no afirmar diferencia de eficiencia menor al 30 %.
- Precedentes directos: Vermeire, Witherden y Vincent (2017), *JCP* 334:497–521 ****, PyFR
  contra STAR-CCM+ igualando **grados de libertad**, no mallas; Coppeans, Fidkowski y Martins
  (2023), AIAA 2023-1845 ****; Jia et al. (2019), *AIAA J.* 57:1636–1648 ****.
- El precedente más cercano a esta tesis: Jiang y Cheng (2021), *Phys. Fluids* 33 ****,
  Nektar++ con SVV contra OpenFOAM v5.0 sobre cilindro; el propio proyecto Nektar++ reporta
  ~30 000 horas de núcleo contra ~1 000 000 estimadas para OpenFOAM al mismo nivel de convergencia
  ****. **Caveat obligatorio: la cifra de OpenFOAM es una extrapolación, no una medición.**
- Qué queda abierto: no hay comparación de costo a error fijo sobre geometría de aplicación en
  régimen supersónico externo.

### 2.7 El hueco identificado

Formularlo como intersección de tres carencias verificadas:

1. El operador SVV no existe en el solver compresible de Nektar++, y la infraestructura que lo
   haría barato existe pero está desconectada.
2. La literatura de SVV compresible valida la formulación pero no la física: todos los casos
   publicados corren a Mach efectivamente incompresible.
3. No hay evidencia de costo a error fijo frente a FVM de segundo orden sobre geometría de
   aplicación en este régimen, y el propio taller de referencia declara que para soluciones no
   suaves *"there is not enough evidence"*.

---

## Capítulo 3 — `OBJETIVOS_Y_ALCANCE`

### 3.1 Pregunta de investigación

**Pregunta principal.** ¿De qué forma se puede lograr, con un solver basado en métodos
espectrales, la misma robustez y precisión que un solver FVM estándar, a menor costo de cómputo,
en flujo compresible externo con M = 0.3 a 3?

**La propuesta.** Un operador de viscosidad espectral desvaneciente modulado por sensor modal,
como mecanismo único de estabilización y robustez.

**Por qué hace falta agregarlo: lo que el `CompressibleFlowSolver` tiene hoy no alcanza.**
Cinco razones, todas verificadas contra el código o contra la literatura primaria. Ésta es la
justificación que sostiene la tesis entera; conviene que esté escrita con este nivel de respaldo.

1. **Sin mecanismo adicional, el esquema no llega al final de la corrida.** Rubio et al. (2025)
   **** reportan divergencia de DG estándar sobre TGV a Re = 1600 en `t/t_c ≈ 4` con flujo
   central y `≈ 7` con Roe: *"The Riemann solver alone does not provide sufficient dissipation to
   stabilize the simulation."* La disipación de upwind, que es el mecanismo por defecto del
   solver, es insuficiente en el régimen sub-resuelto.
2. **Los mecanismos existentes son dos dispositivos separados con constantes empíricas, y sus
   propios autores lo dicen.** Persson y Peraire (2006) ****: `κ` se elige empíricamente; a
   Mach alto hay que subir el coeficiente y eso ensancha el choque; el sensor *"is not very
   effective at discriminating between shocks and contact discontinuities"*; el coeficiente *"has
   been tuned for each Mach number"*. Moxey et al. (2020) describen el `μ_av` de Nektar++ como
   proporcional a un *"O(1) constant"* sin valor derivado **** — y publican `μ0 = 0.25`
   (Yan et al., 2021 ****) mientras el código trae `1.0` por defecto ****.
3. **La viscosidad artificial existente es constante por elemento: puede fijar cuánta disipación
   introduce, pero no dónde en el espectro.** Ése es precisamente el grado de libertad que separa
   estabilizar de arruinar el resultado en régimen sub-resuelto, y es lo que Rubio et al. (2025)
   identifican como carencia: *"excessive SGS dissipation at intermediate scales can corrupt the
   solution near the cutoff wavenumber, especially in high-order discretizations where numerical
   dissipation is minimal [ These findings emphasize the need for adaptive stabilization
   strategies both in space and in time. An example is spectral vanishing viscosity"* ****.
   Es una recomendación explícita, por autores del campo, de hacer exactamente esto.
4. **Los mecanismos existentes tienen huecos de cobertura verificados en el código** ****:
   la viscosidad artificial física no está disponible para Euler; el sensor de dilatación y el
   filtro de Ducros sólo son alcanzables con `ShockCaptureType = Physical`; **el filtro modal
   exponencial está implementado sólo para segmento, cuadrilátero y hexaedro** — sobre triángulos,
   tetraedros, prismas o pirámides aborta con `ASSERTL0`. Para la malla del Aconcagua, que será
   híbrida con prismas de capa límite, ese mecanismo directamente no corre.
5. **Ninguno de los mecanismos existentes tiene teoría de convergencia; SVV sí.** El resultado de
   Tadmor (1989) **** garantiza convergencia a la solución de entropía sin pérdida de precisión
   espectral — para leyes escalares 1D, y hay que decir ese alcance. Ningún resultado análogo
   existe para viscosidad artificial con sensor.

**Y la razón de oportunidad**, que no es teórica sino de ingeniería: los núcleos SVV ya están
implementados en `library/StdRegions` para **todas** las formas de elemento, con tres variantes
de núcleo (exponencial, ley de potencias, núcleo DG de Moura) ****. El costo marginal de
llevarlos al solver compresible es bajo comparado con el valor. La brecha es de conexión, no de
desarrollo.

**Sub-preguntas**, cada una mapeada a una sección de resultados:

- ¿Qué combinación de amplitud y corte estabiliza sin dominar la disipación física?
- ¿Cómo depende esa combinación de `M`? ¿Es una ley suave o hay transición de régimen?
- ¿Es el resultado igual o mejor que el de los mecanismos ya disponibles en el solver, al mismo
  costo?
- ¿Cuál es la curva costo–error frente a OpenFOAM sobre el Aconcagua?

### 3.2 Objetivo general

Implementar y calibrar un operador SVV modulado por sensor modal en el `CompressibleFlowSolver`
de Nektar++, de modo que provea estabilización única desde régimen subsónico turbulento
sub-resuelto hasta supersónico con choques en M = 0.3–3, y medir el costo a precisión fija del
esquema resultante frente a un solver FVM de segundo orden sobre un caso de aplicación con
referencia experimental.

### 3.3 Objetivos específicos

1. Inventariar los mecanismos de estabilización del solver compresible, sus parámetros y sus
   límites documentados, y determinar el punto de inserción del operador. *(Hecho en parte; ver
   cap. 8.)*
2. Reconstruir y documentar la arquitectura del solver compresible hasta el nivel de detalle que
   permita insertar el operador, y exponerla como resultado del trabajo: el mapa de cómo se ensambla
   el término difusivo en Nektar++ es en sí mismo material publicable, porque no existe escrito en
   ningún lado con ese nivel de detalle.
3. Implementar el operador como subclase de `ArtificialDiffusion`, registrada en su factory.
4. **Establecer un tratamiento de trazas consistente con la formulación DG.** *Qué es esto:* el
   término SVV es difusivo, de segundo orden, y en DG un término de segundo orden no se puede
   discretizar sin decidir qué valor toma la solución —y su gradiente— en las caras donde los
   elementos no coinciden. La receta habitual introduce una variable auxiliar `q = ∇u` y define
   dos flujos numéricos de interfaz: uno para `q` y otro para el flujo viscoso. Eso es el
   tratamiento de trazas. La decisión concreta acá es **dónde entra el filtro modal en esa
   cadena**: se puede filtrar el gradiente `q` después de que el operador LDG lo reconstruyó
   —que es lo que permite la arquitectura existente, porque `ArtificialDiffusion` ya trae su
   propio operador LDG y entrega `qfield`— o filtrar antes y redefinir el flujo de interfaz. La
   primera opción reutiliza maquinaria ya verificada; la segunda es más limpia formalmente pero
   exige reescribir flujos. **El riesgo asociado es menor de lo que parecía en v0.2**, justamente
   porque la clase base ya resuelve la parte difícil: la modificación es una operación modal por
   elemento sobre `qfield`.
5. Calibrar amplitud y frecuencia de corte en función del número de Mach.
6. Cuantificar la disipación introducida por el operador y compararla con la disipación física,
   usando la descomposición `ε₁` (solenoidal) y `ε₃` (dilatación de presión) del caso TGV.
7. Comparar contra los mecanismos ya disponibles en el solver, al mismo costo y sobre los mismos
   casos.
8. Construir una malla de alto orden del Aconcagua a partir del STL disponible, con curvado de
   superficie y verificación de validez de los elementos.
9. Definir y ejecutar un protocolo de comparación de costo a error fijo contra OpenFOAM, con la
   telemetría de vuelo como referencia experimental.

### 3.4 Criterios de trabajo previos a la comparación

*(Sección corta y de alto valor. Su función es fijar por escrito, antes de correr nada, qué
constituye un resultado aceptable. No es falsación popperiana: es control de la tentación de
ajustar el criterio después de ver el número.)*

- **Criterio de precisión conservada.** Con el operador activo, el orden de convergencia observado
  sobre solución suave (vórtice isentrópico) no debe caer por debajo del orden de diseño menos
  media unidad. Si cae, el operador no es admisible con esos parámetros.
- **Criterio de estabilización.** Existe al menos una combinación `(ε₀, m_N)` con la que Shu–Osher
  y el escalón supersónico M = 3 corren hasta el tiempo final sin oscilaciones espurias visibles
  y sin fallo de positividad.
- **Criterio de balance energético.** En el extremo subsónico, la disipación atribuible al
  operador no debe exceder una fracción declarada de la disipación física total. **Fijar el umbral
  y el estimador acá, antes de correr.** El estimador propuesto es el reportado en la
  especificación del caso TGV del taller (Hillewaert, caso C3.5 ****), que sugiere comparar la
  tasa de disipación de energía cinética consistente contra la variante numérica que incluye los
  términos de salto del esquema DG.
- **Criterio de no inferioridad.** El operador debe igualar o mejorar el resultado del mecanismo
  `NonSmooth` existente, al mismo costo en work units, sobre el conjunto mínimo de casos de 9.3.
  Si sólo iguala, el aporte sigue en pie —por generalidad de forma de elemento y por teoría de
  convergencia—, pero hay que decirlo así.
- **Criterio de comparación de costo.** La comparación contra OpenFOAM se reporta como error
  frente a work units, medidas según la definición de HiOCFD (TauBench, tiempo de solver
  únicamente) ****, y **no se afirma diferencia de eficiencia menor al 30 %** (Wang et al.,
  2013 ****).
- Declarar métrica de error, solución de referencia y rango explorado antes de reportar cualquier
  número.

### 3.5 Alcance y delimitación

**Se incluye:**

- Base de expansión modal jerárquica, la nativa de Nektar++.
- Navier–Stokes compresible como objetivo de implementación.
- Rango M = 0.3 a 3.
- Casos 1D, 2D y 3D con prioridad decreciente por costo, más HB-2 y el Aconcagua.

**Se excluye explícitamente:**

- Desarrollo de formulación matemática nueva.
- **Formulación en variables de entropía y demostración de estabilidad entrópica del esquema.**
  El operador se escribe sobre las **variables conservativas**, que son las que entrega la
  maquinaria de difusión artificial de Nektar++, y se discretiza con el operador LDG que esa
  maquinaria ya provee. No se reformula el solver. Razón verificada, no de conveniencia: tanto
  Mateo-Gabín et al. (2022) como Lin, Chan y Tomas (2023) requieren discretización de volumen por
  *flux differencing* con flujos de dos puntos entrópicamente conservativos y operadores con
  propiedad de suma por partes sobre nodos de Gauss–Lobatto, y Nektar++ no tiene ninguno de esos
  ingredientes ****. **El papel de esa literatura en esta tesis es de justificación**: es lo que
  permite argumentar que el filtro es un mecanismo disipativo válido en régimen compresible, y no
  una heurística. Detalle en 7.3, 7.7 y el Anexo A.1.
- **Preservación garantizada de positividad.** Nektar++ no tiene ningún limitador de positividad
  ****. Los fallos se detectan y se documentan. Un limitador de escalado tipo Zhang y Shu
  (2010) queda como objetivo opcional, no comprometido.
- Régimen hipersónico y efectos de gas real.
- Optimización de la implementación. El costo se mide sobre una implementación correcta y no
  optimizada; declararlo es lo que hace honesta la comparación del capítulo 9.
- LES resuelta en pared del cohete a número de Reynolds de vuelo. Es inviable y está documentado
  que lo es (Slotnick et al., 2014 ****); ver 5.6.

### 3.6 Referencia experimental: telemetría de vuelo del Aconcagua

- Se dispone de telemetría de vuelo del vehículo. Eso convierte el caso de aplicación de
  demostración en validación, que es un salto de categoría — pero sólo si se trata con el método
  correcto.
- **Se usa todo el vuelo, con las suposiciones declaradas.** El tramo de coasteo y el propulsado
  no tienen el mismo estatus y hay que decirlo: en coasteo el empuje es nulo y la masa constante,
  de modo que `C_D` queda determinado por la ecuación de movimiento; en el tramo propulsado hay dos
  incógnitas —empuje y arrastre— y una sola ecuación, así que hace falta una **curva de empuje
  supuesta**, y ésa pasa a ser la suposición dominante. Bawa y Michalski (2025), AIAA 2025-104785
  ****, lo enuncian así: *"In the absence of direct thrust measurements, the reconstruction
  admits an infinite set of pairs, so either thrust or drag must be known for determinacy."*
  Reportan también que el coasteo es ~96 % del ascenso hasta el apogeo, que la reconstrucción es
  notablemente más limpia ahí, y documentan un caso en que un error de la curva de empuje cerca de
  Mach 1.4 produjo un arrastre reconstruido **negativo**.
- **Consecuencia de procedimiento:** las barras de incertidumbre del tramo propulsado son más
  anchas que las del coasteo y hay que graficarlas distinto. Si una discrepancia con el CFD aparece
  sólo en el tramo propulsado, la primera hipótesis es la curva de empuje, no el solver.
- Precedente directo de comparación CFD contra telemetría de cohete universitario: Kierulf (2022),
  ICAS 2022-0414 ****, cuatro vuelos instrumentados; acuerdo entre CFD y `C_D` de vuelo de
  0.7 % en un vuelo y 17 % en otro. **La literatura es delgada y hay que decirlo** — lo que
  también significa que una comparación cuidadosa y con incertidumbre cuantificada es un aporte
  real.
- **Consecuencia de modelado, que con la decisión de usar todo el vuelo pesa más:** un CFD sin
  pluma representa la configuración de coasteo, no la propulsada. Bawa y Michalski documentan
  **** que *"the cessation of thrust combined with the presence of an exposed nozzle produces a
  transient increase in drag, primarily attributed to an increase in base drag."* Es decir: en los
  puntos propulsados, la geometría simulada y la volada difieren en la base. O se modela la pluma
  —lo que es un trabajo en sí mismo y no está en el alcance—, o se declara que la comparación en
  esos puntos tiene un sesgo conocido de arrastre de base, y de qué signo.

### 3.7 Amenazas a la validez

- **Comparación entre un código que se está escribiendo y uno maduro.** Sesgo estructural en
  contra del primero. Mitigación: reportar costo y grados de libertad por separado, perfilar el
  costo del operador SVV aparte del resto del solver, y aplicar la regla del 30 %.
- **Elección de la métrica de error del caso Aconcagua.** Si no hay solución de referencia
  independiente, la referencia es una corrida fina; decir cuál y por qué. El escalón intermedio
  HB-2 (ver 9.7) existe justamente para saber si una discrepancia es del solver o de la telemetría.
- **Incertidumbre de la telemetría.** Fuentes documentadas en las dos referencias de 3.6 ****:
  propiedades másicas (incluyendo que el grano de humo sigue quemando durante el coasteo), ángulo
  de ataque desconocido y fuerzas fuera de eje despreciadas, alineación del IMU con los ejes del
  vehículo, modelo de atmósfera y densidad, alineación temporal entre registradores independientes,
  y ruido y deriva del acelerómetro.
- **Ajuste de parámetros posterior a ver los resultados.** Se previene con 3.4.
- **Transplante de constantes desde la literatura.** Los valores óptimos publicados (`P_SVV = 0.1`,
  `C_S = 0.2`) se obtuvieron con DGSEM de forma split y solver de Riemann de baja disipación
  (Manzanero et al., 2020; Mateo-Gabín et al., 2022 ****). Esta tesis usa DG estándar. Esas
  constantes **no** son transplantables sin recalibración, y decirlo es parte del argumento.

---

# PARTE 2 — TEORÍA

## Capítulo 4 — `MARCO_TEORICO`

**Función:** instalar el aparato que usan los capítulos 5 a 7. Regla de poda: si una ecuación no
reaparece, no va.

### 4.1 Ecuaciones de Navier–Stokes compresibles

- Forma conservativa: continuidad, cantidad de movimiento y energía total; vector de estado y
  tensores de flujo convectivo y difusivo.
- Cierre: gas ideal, ley de Sutherland, ley de Fourier, número de Prandtl.
- Adimensionalización y aparición de `M`, `Re`, `Pr`, `γ`. Declarar la escala — de esto depende la
  interpretación de todos los parámetros del operador.
- Estructura hiperbólico-parabólica y velocidades características; ausencia de proyección de
  presión y su consecuencia sobre la rigidez del sistema a bajo `M`.
- **Entropía matemática y desigualdad de entropía.** Se define acá porque es el criterio que
  selecciona la solución físicamente admisible, y por lo tanto lo que le da sentido al resultado de
  convergencia de Tadmor que fundamenta el operador. Convexidad de la entropía y simetrización del sistema. **No se introduce la formulación en variables de entropía**: el operador de esta tesis actúa sobre las variables conservativas (ver 7.3). Las variables de entropía aparecen sólo al citar la literatura de la que se toma la justificación (7.7) y, si se implementa el diagnóstico opcional, al medir el balance de entropía discreta.

### 4.2 Discretización de Galerkin discontinuo espectral

- Formulación débil elemento a elemento; espacio de aproximación discontinuo; elemento de
  referencia y mapeo geométrico.
- Base modal jerárquica de Karniadakis–Sherwin: modos de vértice, arista y burbuja. **La
  jerarquía es lo que permite hablar de "modos altos" sin ambigüedad**, y por lo tanto lo que hace posible el operador.
- Matriz de Vandermonde y transformación nodal ↔ modal; base ortogonal (`eOrtho_A/B` en Nektar++).
  Es la operación sobre la que se montan tanto el sensor como el filtro. Nota práctica: el problema
  de coordenadas colapsadas en triángulos ya está resuelto en la biblioteca ****.
- Cuadratura de Gauss–Lobatto–Legendre; colocación.
- **Esquemas de forma split y propiedad de suma por partes: qué son y por qué esta tesis no los
  tiene.**
  Sección corta pero necesaria, porque es la que justifica la exclusión de 3.5 y prepara el
  Anexo A.1. Referencias: Gassner, Winters y Kopriva (2016), *JCP* 327:39–66 ****;
  Chan (2018), *JCP* 362:346–374 ****.
- Flujos numéricos de interfaz: Roe, HLLC, Lax–Friedrichs local, AUSM. Rol del upwinding como
  fuente de disipación *ya presente* — hay que contabilizarla al medir la del operador. Nota
  verificada ****: ninguno de los flujos disponibles en Nektar++ es entrópicamente estable por
  construcción.
- Tratamiento del término difusivo: LDG y penalización interior (IP). Ésta es la maquinaria que
  el capítulo 7 reutiliza para las trazas del término SVV.
- Integración temporal explícita Runge–Kutta, restricción CFL y su escalamiento con `P²`. Nota:
  el término SVV es difusivo y endurece esa restricción; cuantificarlo es parte del costo.

### 4.3 Propiedades del esquema

- Convergencia espectral en soluciones suaves; degradación a orden algebraico ante discontinuidad.
- **Análisis de dispersión y disipación.** Es el resultado que hace que la elección del corte no
  sea arbitraria. Moura, Sherwin y Peiró (2015), *JCP* 298:695–710 **** proponen un criterio
  para estimar la resolución efectiva de DG en términos del mayor número de onda que resuelve con
  precisión, y lo validan sobre turbulencia de Burgers. Moura et al. (2017), *JCP* 330:615–623
  **** lo extienden a espectros 3D como la "regla del 1 %". Alcance honesto: la teoría es
  lineal, la validación es a posteriori.
- Análisis no modal: Fernandez, Moura, Mengaldo y Peraire (2019), *CMAME* 346:43–62 ****,
  que plantea la pregunta en los términos exactos de esta tesis: *"why do SEM sometimes suffer
  from numerical stability issues in LES? And, why do they at other times be robust and
  successfully predict under-resolved turbulent flows even without a subgrid-scale model?"*
- Aliasing: origen en la cuadratura de productos no lineales, desaliasing por sobre-integración,
  y su distinción de la estabilización por SVV. **No son lo mismo y confundirlos es un error
  frecuente.** Kirby y Sherwin (2006) advierten explícitamente **** que SVV no repara los
  errores de integración inconsistente. Además, la propia documentación de Nektar++ reconoce
  **** que *"we are using a discontinuous discretisation, meaning that aliasing effect are not
  fully controlled, since the boundary terms introduce non-polynomial functions into the problem."*

---

## Capítulo 5 — `TURBULENCIA_Y_ESCALAS`

**Función:** decir qué escalas se resuelven, cuáles no, qué se supone que pasa con las que no se
resuelven, y en qué régimen físico está cada caso. Sin este capítulo, la palabra
"disipación" queda ambigua durante toda la tesis: no se distingue la física de la numérica.

### 5.1 Escalas de la turbulencia y qué significa resolverlas

- Cascada de energía, escalas integral y de Kolmogorov, y el conteo de grados de libertad que
  hace inviable la DNS a `Re` de vuelo.
- Definición operativa de "resuelto": el criterio de resolución efectiva de Moura et al. (2015,
  2017) **** — el número de onda hasta el cual el esquema transporta sin disipación numérica
  apreciable. Es la definición que esta tesis usa, porque es la que se puede medir.
- Número de onda de Nyquist del esquema, `k_Ny = π/Δx` con `Δx = L/[(P+1)N_el]` ****.

### 5.2 LES compresible: qué cambia respecto del caso incompresible

- Filtrado de Favre `f̃ = ρf/ρ̄`: mantiene la continuidad en su forma original.
- Tensor de esfuerzos submalla en cantidad de movimiento.
- **La ecuación de energía es donde proliferan los términos**: flujo de calor submalla, difusión
  turbulenta submalla, difusión viscosa submalla, dilatación de presión submalla y disipación
  submalla. No hay agrupamiento canónico único; distintos autores parten distinto. Referencias:
  Erlebacher et al. (1992) ****; Moin et al. (1991), *Phys. Fluids A* 3:2746–2757
  ****; Vreman, Geurts y Kuerten (1995), *Appl. Sci. Res.* 54:191–203 **** para la
  magnitud relativa de cada término; Garnier, Adams y Sagaut (2009), cap. 3 **** como
  revisión.
- **Advertencia:** ninguno de esos textos fue leído para este documento. No atribuirles rankings
  de términos ni números sin abrir la fuente.

### 5.3 Compresibilidad de la turbulencia: número de Mach turbulento

- Definición `M_t = √(u'ᵢu'ᵢ)/c̄`.
- Disipación dilatacional: Sarkar et al. (1991) **** proponen `ε_d = α₁ M_t² ε_s`; Zeman (1990)
  propone una alternativa. Ambos calibrados sobre **turbulencia homogénea y capas de mezcla
  compresibles**.
- **El resultado que decide el alcance de esta tesis.** Duan, Beekman y Martín (2011), *JFM*
  672:245–267 ****: DNS de capas límite turbulentas con Mach de corriente libre de 0.3 a 12,
  *"the explicit dilatation terms such as pressure dilatation and dilatational dissipation remain
  small for the present Mach number range"*, y la transformación de van Driest, la relación de
  Walz, el escalado de Morkovin y la analogía fuerte de Reynolds siguen valiendo.
- **Enunciado honesto para el cohete:** en capa límite adherida sobre cuerpo esbelto a M = 1.5–3,
  la hipótesis de Morkovin vale; densidad y viscosidad varían fuertemente, pero la estructura de
  la turbulencia no. Las correcciones de Sarkar y Zeman **no** hacen falta. La compresibilidad que
  importa es la del flujo medio: densidad variable, choques, abanicos de expansión, interacción
  choque–capa límite y flujo de base.
- **Caveat que hay que escribir en el mismo párrafo:** eso vale para capas *adheridas*. En la
  estela de base del cohete se forma una capa de corte supersónica libre, y ahí `M_t` y el número
  de Mach convectivo vuelven a ser relevantes. Si el caso incluye la base, declarar la división.

### 5.4 LES implícito en DG: "el esquema numérico es el modelo"

- La afirmación fuerte: Beck et al. (2014), *IJNMF* 76:522–548 ****, *"Without the need for any
  additional filtering, explicit or implicit modelling, or artificial dissipation, our high-order
  schemes capture the turbulent flow at the considered Reynolds number range very well"* — sobre
  cilindro a Re = 3900, colina periódica y SD7003. Atribuyen el mecanismo al desaliasing polinomial
  y al orden alto, no a un sustituto de viscosidad turbulenta.
- Su fundamento y el límite de ese fundamento: análisis lineal de dispersión y disipación
  (Moura et al., 2015, 2017 ****), extrapolado por analogía a turbulencia 3D no lineal y
  validado a posteriori sobre Burgers y TGV. **No es una derivación de un cierre submalla, y no
  afirma fidelidad física de la transferencia submalla** — sólo que el rango de disipación
  numérica se puede ubicar donde termina el rango resuelto.
- Gassner y Beck (2013), *TCFD* 27:221–237 ****: cita universal para esta capacidad;
  texto completo no accesible, no citar afirmaciones puntuales sin leerlo.
- El matiz de orden: Winters et al. (2018), *JCP* 372:1–21 ****, mejor calidad a orden
  moderadamente alto (~6), y órdenes muy altos muestran *"spurious features attributed to a
  sharper dissipation in wavenumber space"*.

### 5.5 Dónde se ubica SVV entre LES implícito y LES explícito

- SVV como modelo LES: la idea es de Karamanos y Karniadakis (2000) **[V-bib, texto no leído]**,
  según lo caracterizan Pasquetti (2006) y Kirby y Sherwin (2006), ambos leídos ****. **Toda
  esa línea es incompresible.**
- **No existe equivalencia demostrada con un modelo de viscosidad turbulenta.** Pasquetti (2006)
  **** es explícito: los escalados teóricos de los parámetros están establecidos para leyes
  escalares, *"but up to our knowledge such results are not yet available for the Navier-Stokes
  equations"*. Lo más cercano a una equivalencia es de *disipación numérica*, no de modelo físico:
  Moura et al. (2016) **** muestran que CG con SVV tiene disipación *"broadly equivalent"* a la
  del upwinding en DG.
- El vínculo con Smagorinsky es **constructivo, no derivacional**: Manzanero et al. (2020) ****
  fijan la amplitud con la fórmula de Smagorinsky y después la filtran espectralmente. Elimina un
  parámetro libre; no deriva una equivalencia.
- El resultado físico que hace atractivo el híbrido, de Manzanero et al. (2020) ****: *"the main
  difference with the standard Smagorinsky model is that the latter presents an excessive
  dissipation when the flow is laminar [ This is naturally avoided with the SVV technique, as it
  filters-out the laminar (smooth) energy components of the dissipation."*
- **La síntesis defendible, y es la frase que hay que escribir:** la disipación numérica de DG se
  puede caracterizar lo bastante bien como para ubicar un rango de disipación artificial en el
  límite de resolución, y en LES bien resuelta eso alcanza sin modelo submalla explícito. No es un
  cierre físico, depende de la resolución y del esquema, y en régimen severamente sub-resuelto es
  insuficiente (DG estándar) o está mal ubicado en número de onda. SVV es el mecanismo que permite
  fijar la **cantidad** y la **ubicación espectral** de esa disipación de forma independiente — y
  por eso se ubica entre LES implícito y LES explícito, sin ser ninguno de los dos.
- Respaldo contemporáneo de esa lectura: Rubio et al. (2025) ****, ya citado en 3.1.

### 5.6 Estrategia de resolución de escalas para cada caso

- **TGV a Re = 1600:** es formalmente una DNS (especificación C3.5 del taller ****); corrida
  sub-resuelta se convierte en LES implícito o DNS sub-resuelta. Decir cuál se está haciendo.
- **Casos 1D y 2D con choques:** no hay turbulencia; el operador se evalúa como mecanismo de
  captura, no como modelo.
- **Cohete Aconcagua:** acá hay que ser explícito y no prometer de más.
  - Para capas límite adheridas a `Re` de vuelo nadie corre LES. La estimación de Spalart citada
    en Slotnick et al. (2014) **** ubica la LES de ala completa para uso ingenieril en 2045.
  - En la literatura de lanzadores, la resolución de escalas se usa **selectivamente en las zonas
    separadas** —base, estela de retro-cuerpo, interacción choque–capa límite—, con RANS/LES zonal
    o DES: Statnikov et al. (2015), *Phys. Fluids* 27:016103 ****; Simon et al. (2006),
    *AIAA J.* 44:2578–2590 ****; Deck (2012), *TCFD* 26:523–550 ****.
  - **Lo que el caso del cohete puede afirmar honestamente:** que el esquema estabilizado por SVV
    se mantiene estable y produce resultados físicamente sensatos en configuración externa
    supersónica con choques y turbulencia sub-resuelta; y que la curva costo–error es la que es.
    **Lo que no puede afirmar:** capa límite turbulenta resuelta, LES validada, o cargas de buffet
    predictivas.
  - Un párrafo sobre la trampa de la transición: Slotnick et al. (2014) **** reportan que la
    región laminar y transicional requiere de 10 a 100 veces más celdas que la turbulenta, lo que
    hace que un cuerpo esbelto de proa afilada sea *más* difícil, no más fácil, que la intuición.

---

## Capítulo 6 — `FENOMENOS`

**Función:** justificar el operador por la vía del problema. Cada sección describe un modo de
falla, lo caracteriza espectralmente y enuncia qué le exige al mecanismo estabilizador. La
conclusión es que las exigencias son compatibles entre sí, y por eso un solo operador puede
cubrirlas.

### 6.1 Sub-resolución y acumulación de energía en modos altos

- Manifestación discreta: acumulación de energía en el extremo alto del espectro modal, porque el
  esquema no tiene mecanismo de transferencia hacia escalas que no existen.
- Consecuencia verificada, no supuesta: Rubio et al. (2025) **** documentan la divergencia de
  DG estándar sobre TGV Re = 1600 entre `t/t_c ≈ 4` y `≈ 7`.
- Exigencia al mecanismo: disipación selectiva en el extremo alto, nula en los modos resueltos,
  de magnitud comparable a la transferencia submalla física.

### 6.2 Discontinuidades y fenómeno de Gibbs

- Tipos presentes en el rango: choque normal y oblicuo, discontinuidad de contacto, abanico de
  expansión.
- Oscilaciones de Gibbs: amplitud que no decae con el refinamiento, contaminación dentro del
  elemento. Referencia: Gottlieb y Shu **[F]**, verificar antes de citar.
- Consecuencias catastróficas: presión o densidad negativas, fallo de la ecuación de estado.
- Exigencia: disipación **localizada** y de magnitud suficiente — es decir, activación por sensor,
  no uniforme. Evidencia de que SVV solo no alcanza en choque fuerte: Mateo-Gabín et al. (2022)
  ****, *"SVV filtering is not enough to control the oscillations near the strong shock."*
  **Esto justifica el sensor como parte no opcional de la propuesta.**

### 6.3 Inestabilidad por aliasing

- Cómo un error de cuadratura en los flujos no lineales realimenta energía a modos altos.
- Por qué se confunde con 6.1 y cómo distinguirlos en la práctica: respuesta al desaliasing por
  sobre-integración frente a respuesta al filtro.
- Límite estructural en DG, documentado por Nektar++ **** (citado en 4.3).

### 6.4 Interacción choque–turbulencia

- Régimen M ≳ 1.5 con turbulencia incidente; amplificación de fluctuaciones al atravesar el
  choque.
- Por qué es el caso más exigente: los dos modos de falla anteriores ocurren en el mismo elemento,
  y un esquema con dos mecanismos separados debe decidir cuál aplica. **Es el argumento más fuerte
  a favor del operador único: elimina la decisión.**

### 6.5 Pérdida de positividad

- Condición de admisibilidad del estado; dónde falla típicamente: base del cuerpo, expansión
  fuerte, arranque impulsivo.
- Estado en el código: Nektar++ no tiene ningún limitador de positividad ****. Se declara como
  fenómeno observado y documentado, no corregido.

### 6.6 Síntesis: requisitos del mecanismo estabilizador

- Tabla de requisitos derivados: selectividad espectral, localidad espacial, activación por
  suavidad, consistencia (se anula bajo refinamiento), preservación del orden en zonas suaves,
  costo marginal bajo, cobertura de todas las formas de elemento, y compatibilidad con la
  formulación DG existente.
- Contrastar contra los mecanismos disponibles del capítulo 8 y mostrar cuáles satisfacen qué.
  **El filtro exponencial cae por cobertura de forma de elemento** ****; la viscosidad
  artificial cae por selectividad espectral; el desaliasing no es un estabilizador. Ésta es la
  bisagra del capítulo 7.

---

## Capítulo 7 — `FILTRO_SVV`

**Reestructurado en v0.3.** El orden ahora es: qué garantiza la teoría → qué formulación se
adopta → qué de esa formulación es alcanzable en este código y qué no → qué se mide en lugar de
demostrar.

### 7.1 Qué garantiza la teoría, y hasta dónde

- Por qué la viscosidad artificial clásica destruye la precisión: opera sobre todos los modos.
- Idea de Tadmor: operador espectralmente selectivo cuya amplitud se desvanece con la resolución.
- **Enunciado preciso del resultado y de su alcance.** Tadmor (1989) ****: convergencia a la
  solución de entropía sin pérdida de precisión espectral, para **leyes de conservación escalares
  no lineales en una dimensión, base de Fourier**. Maday, Ould Kaber y Tadmor (1993) ****
  extienden a base de Legendre y a dinámica de gases 1D. **No hay resultado análogo para sistemas
  no lineales multidimensionales.** Decirlo es obligatorio y es barato.

### 7.2 Formulación del operador

- Forma general: término adicional `∂ₓ(ε_N Q_N ∂ₓu)`, con `Q_N` diagonal en la base modal
  ortogonal.
- Núcleo de Fourier frente a núcleo de Legendre, y por qué el segundo es el aplicable a elementos.
- Forma spectral/hp de Kirby y Sherwin (2006) ****: núcleo exponencial, frecuencia de corte
  `m_N`, amplitud `ε_N`, escalamiento con `h` y `P`.
- **Tres variantes de núcleo, todas ya implementadas en la biblioteca** ****: exponencial
  (`eFactorSVVCutoffRatio` + `eFactorSVVDiffCoeff`), ley de potencias de Moura
  (`eFactorSVVPowerKerDiffCoeff`), y núcleo DG de Moura (`eFactorSVVDGKerDiffCoeff`). La
  comparación entre las tres es casi gratuita y es material de tesis, no sobrecarga.
- Argumento a favor del núcleo DG: Moxey et al. (2020) **** lo describen como uno que
  *"replicates the desirable dispersion and diffusion properties of DG schemes and does not
  require the manual tuning of parameters found in the classical SVV formulation."*
- Interpretación: filtro modal aplicado a la derivada, no a la solución. La distinción importa
  para la conservación.

### 7.3 Extensión al sistema compresible

- **Decisión tomada: el operador actúa sobre las variables conservativas.** Es la formulación
  nativa de Nektar++ y la que permite entrar sin reformular nada: `ArtificialDiffusion` ya
  construye su operador LDG sobre las variables conservativas y entrega el gradiente `qfield`
  ****, de modo que el operador se reduce a una operación modal por elemento sobre ese
  gradiente. Es también la decisión que mantiene el trabajo dentro de las 16 semanas disponibles.
- **Qué se pierde con esa elección, dicho explícitamente.** La demostración de disipación no
  negativa del filtro de Mateo-Gabín et al. (2022) —Teorema 1 y Propiedad 4, más la factorización
  de Cholesky de su sec. 5.2— está enunciada **en variables de entropía** ****. Sobre variables
  conservativas no se transfiere literalmente. Esta tesis no puede reclamarla como propiedad
  demostrada de su implementación; la usa como argumento de plausibilidad, y lo dice así (ver 7.7).
- **Acoplamiento mediado por la ecuación de estado, que con esta elección importa más, no menos.**
  Filtrar densidad y energía total por separado puede producir un estado termodinámicamente
  inconsistente —presión o temperatura sin sentido físico— si no se cuida la combinación. Es el
  punto donde la extensión compresible deja de ser trivial y hay que desarrollarlo, no mencionarlo.
- **Consecuencia práctica a verificar en los primeros casos 1D:** si el filtrado por componentes conservativas
  genera estados inadmisibles en las zonas de gradiente fuerte, la alternativa dentro del mismo
  alcance es filtrar sobre variables primitivas o sobre un subconjunto de componentes. Dejarlo
  planteado como bifurcación, con el criterio de decisión escrito de antemano.
- Forma matricial del término añadido y su lugar dentro del tensor de flujo viscoso, que es lo que
  permite reutilizar la maquinaria LDG existente.
- Escalamiento de la amplitud con el número de Mach: motivación desde el capítulo 5, forma
  funcional propuesta `ε_N = ε₀ f(M)`. **El mapa paramétrico del objetivo 5 es la determinación
  empírica de `f`.** Nota de honestidad: la literatura no ofrece `f`, porque toda ella corre a
  Mach ≈ 0.1.

### 7.4 Sensor modal de suavidad

- Definición de Persson y Peraire (2006) ****: `s_e = log₁₀(‖q − q̃‖²/‖q‖²)`, con `q̃` la
  proyección de orden reducido, y rampa senoidal con umbral `s₀ = S_κ − 4.25 log₁₀(p)`.
- Implementación en Nektar++ ****: `VariableConverter::GetSensor`; parámetros `Skappa`
  (por defecto −1.0), `Kappa` (0.25), `SensorOffset` (1). Nota verificada: devuelve cero para
  elementos con `numModes ≤ SensorOffset`, de modo que elementos de orden 0 y 1 no reciben nada.
- Limitaciones documentadas por sus autores ****, ya citadas en 3.1: umbral empírico, mala
  discriminación entre choque y contacto, necesidad de subir el coeficiente a Mach alto.
- Sobre qué variable se evalúa (densidad, presión, o indicador combinado) y por qué la elección
  cambia el comportamiento en contactos frente a choques.
- **Papel dual del sensor en esta tesis**, que es la forma concreta que toma "un solo mecanismo":
  modula el operador entre el régimen de amplitud baja siempre activa (estabilización de
  sub-resolución) y el de amplitud alta localizada (captura de choque). Explicitarlo.

### 7.5 Selección de parámetros

- Criterios analíticos de Moura, Sherwin y Peiró (2016) ****: qué predice el análisis de
  autosoluciones sobre `(ε_N, m_N)`, y la corrección del número de Péclet local.
- **Predicción previa frente a calibración empírica.** Declarar el rango esperado antes del
  barrido y luego contrastar. Eso convierte el capítulo 10 en verificación y no en ajuste.
- Valores de referencia de la biblioteca, como punto de partida documentado ****:
  `SVVCutoffRatio = 0.75` (el primer 75 % de las frecuencias no se amortigua), `SVVDiffCoeff = 0.1`
  para núcleo exponencial y de potencias, `1.0` para el núcleo DG.
- **Detalle de diseño que hay que discutir explícitamente** ****: en el solver incompresible el
  coeficiente SVV se inyecta normalizado por la difusividad física (`m_sVVDiffCoeff / m_kinvis`).
  En un solver compresible no hay una `Kinvis` única y `μ` varía espacialmente, de modo que esa
  normalización no tiene traducción directa. Es un punto de formulación real, no un detalle de
  implementación.
- Valores publicados en DG compresible: `P_SVV ≈ 0.1` y `C_S = 0.2` (Manzanero et al., 2020;
  Mateo-Gabín et al., 2022 ****), con la advertencia de no transplante de 3.7.
- Restricción de consistencia energética y su estimador (ver 3.4).
- Impacto sobre el paso temporal: cuantificarlo, es parte del costo reportado en el capítulo 9.

### 7.6 Tratamiento de trazas en formulación DG

- El problema, explicado como en 3.3.4: el término es difusivo y requiere flujos numéricos de
  interfaz; una elección inconsistente rompe conservación o estabilidad.
- Opciones: filtrar `qfield` después del operador LDG que ya provee la clase base —ruta de menor
  riesgo, porque reutiliza flujos verificados—, o redefinir el flujo de interfaz sobre el gradiente
  filtrado.
- Criterio de decisión, y verificación por vía de los tests unitarios y del caso de orden de
  convergencia: si el par de flujos heredado deja de ser consistente al intercalar el filtro, el
  síntoma es la caída del orden observado en flujo suave, no un fallo ruidoso.
- Nota de riesgo revisada: **menor de lo estimado en v0.2**, porque `ArtificialDiffusion` ya
  construye su propio operador LDG y entrega `qfield` ****; la modificación es una operación
  modal por elemento.

### 7.7 Por qué el filtro es defendible en compresible: el papel de la literatura entrópica

*(Ésta es la sección que sostiene la validez del mecanismo. No se implementa ninguna formulación
entrópica: se la usa como argumento. Desarrollo completo en el Anexo A.1.)*

**El problema que resuelve esta sección.** El operador tiene teoría de convergencia sólo para
leyes de conservación escalares en una dimensión (7.1). Aplicarlo a Navier–Stokes compresible es
una extrapolación. La pregunta legítima del tribunal es: ¿por qué habría de funcionar, más allá de
que empíricamente funcione? La respuesta se construye en cuatro pasos, ninguno de los cuales
requiere reformular el solver.

1. **El mecanismo es disipativo por construcción, no por accidente.** Mateo-Gabín et al. (2022)
   demuestran **** que para cualquier núcleo de filtro positivo el término filtrado tiene
   contribución no negativa a la disipación, vía una factorización de Cholesky de la matriz viscosa.
   El resultado está enunciado en variables de entropía y por lo tanto **no se transfiere
   literalmente** a esta implementación (7.3) — pero muestra que la extensión de SVV al sistema
   compresible admite una demostración bajo una formulación vecina, y eso es cualitativamente
   distinto de una heurística sin respaldo. Enunciarlo con esa precisión, sin sobrevender.
2. **Existe evidencia empírica directa en DG compresible.** Manzanero et al. (2020) y Mateo-Gabín
   et al. (2022) **** hacen funcionar SVV sobre las ecuaciones compresibles, con choques
   (Shu–Osher, escalón M = 3) y en régimen turbulento sub-resuelto (Taylor–Green). Caveat que hay
   que escribir en el mismo párrafo: ambos corren a Mach ≈ 0.1 y sobre esquemas de forma split con
   solver de Riemann de baja disipación, de modo que sus constantes no son transplantables (3.7).
3. **La alternativa —no poner nada— está descartada por evidencia.** Rubio et al. (2025) ****
   documentan que DG estándar sin forma split diverge sobre Taylor–Green Re = 1600 entre
   `t/t_c ≈ 4` y `≈ 7`. Es la configuración de esta tesis.
4. **Lo que reemplaza a la demostración es la medición.** Diagnóstico opcional pero barato: balance
   de entropía discreta `∫ρs` y su derivada temporal, mínimos de densidad y presión por paso.
   Comparar DG puro, DG con la viscosidad artificial existente y DG con SVV sobre los mismos casos.
   Un gráfico de `dS/dt` que se va positivo en el caso base y se mantiene acotado con SVV es un
   resultado propio y honesto, y no exige ninguna reformulación.

**Y la ruta que queda mapeada para después.** Chan (2025), *JCP* 543:114380 ****, prueba que
*"a standard DG method with entropy correction artificial viscosity satisfies the same global
entropy inequality satisfied by flux differencing entropy stable DG methods"*, siempre que los
flujos de interfaz se calculen con la proyección de entropía. Es un término viscoso tipo BR1 —
arquitectónicamente la **misma forma** que el operador de esta tesis. Caveats: sus experimentos son
Euler, no Navier–Stokes; y el propio autor advierte que esa viscosidad *"is not intended to damp
spurious oscillations"*, es decir que complementa a SVV en lugar de reemplazarlo.

---

*Material de respaldo de la sección, para el Anexo A.1: qué exigiría adoptar la formulación
entrópica y por qué está fuera de alcance.*

- **Qué exige la formulación entrópicamente estable y con preservación de positividad.** Lin, Chan
  y Tomas (2023), *JCP* 475:111850 ****: *flux differencing* con flujos de dos puntos
  entrópicamente conservativos (usan el de Chandrashekar), operadores de suma por partes de norma
  diagonal, cuadratura de Gauss–Lobatto, un segundo juego de operadores de bajo orden dispersos,
  mezcla por elemento tipo Zhang–Shu, y paso temporal SSP. Mateo-Gabín et al. (2022) **** exigen
  lo mismo del lado inviscido: DGSEM de Gauss–Lobatto con SBP-SAT y forma split.
- **Qué tiene Nektar++.** Nada de eso, verificado por búsqueda exhaustiva en el árbol de fuentes
  v5.10.0 ****: `AdvectionWeakDG` en forma débil clásica, diez solvers de Riemann ninguno
  entrópicamente estable, sin variables de entropía, sin operadores SBP, sin limitador.
- **Conclusión.** Adoptar la formulación literalmente exige reescribir primero la discretización de
  volumen del solver. Está fuera de alcance, y decirlo con esta evidencia es más fuerte que
  omitirlo: la tesis no evita la formulación entrópica por desconocimiento, sino por una
  restricción del código anfitrión que puede nombrar con precisión.
- **Positividad.** Nektar++ tampoco tiene limitador ****. El de escalado de Zhang y Shu
  (2010, 2011) **** es el único de la familia que no presupone forma split, y queda como
  objetivo opcional, no comprometido.
- **Lo que esto le da a la tesis:** poder decir, con citas y con verificación de código, *por qué*
  la ruta de *flux differencing* está cerrada en Nektar++ y *cuál* formulación publicada la reabre.
  Eso es lo que convierte un Proyecto Final sólido en una propuesta de doctorado.

---

# PARTE 3 — IMPLEMENTACIÓN

## Capítulo 8 — `BIBLIOTECA_NEKTAR`

**Función:** que un lector con el código a mano pueda reproducir la implementación. Escribir contra
versión fija y decir cuál. *(Los hallazgos de esta sección están verificados contra Nektar++
v5.10.0, GitLab oficial, commit `33a2a9a`, agosto 2026.)*

### 8.1 Arquitectura de Nektar++

- Jerarquía: `LibUtilities`, `StdRegions`, `LocalRegions`, `SpatialDomains`, `MultiRegions`,
  `SolverUtils`.
- `CompressibleFlowSolver`: `CompressibleFlowSystem::DoOdeRhs` ensambla
  `−Advection(WeakDG + Riemann) + Diffusion(LDGNS o InteriorPenalty) + Forcing`. Patrón factory de
  registro de componentes.
- Archivo de sesión XML: qué se parametriza y qué se compila.
- Restricciones documentadas ****: sólo `WeakDG` está plenamente soportado como operador
  advectivo; sólo `LDGNS` e `InteriorPenalty` como difusivos; los demás funcionan únicamente sobre
  cuadriláteros.

### 8.2 Inventario de la estabilización existente y sus límites

*Tabla, una fila por mecanismo, con columnas: parámetro de sesión | dónde vive en el código | qué
hace | límite documentado | fuente.* Contenido verificado:

| Mecanismo | Parámetro | Límite verificado |
|---|---|---|
| Viscosidad artificial laplaciana | `ShockCaptureType = NonSmooth` | Sólo solvers explícitos; no admite sensor de Ducros ni de dilatación; no admite suavizado C0; fuerza operador LDG |
| Viscosidad artificial física | `ShockCaptureType = Physical` | Sólo Navier–Stokes, no Euler; los dos sensores tienen escalados distintos en `h/p` sin documentar |
| Sensor modal | `ShockSensorType = Modal` | Constante por elemento; umbrales empíricos; devuelve cero en elementos de orden bajo; la guía dice que es el valor por defecto y el código trae `Dilatation` |
| Sensor de dilatación / Ducros | `ShockSensorType`, `DucrosSensor` | Sólo alcanzables con AV física |
| Filtro modal exponencial | `ExponentialFiltering` | **Implementado sólo para segmento, cuadrilátero y hexaedro**; aborta en triángulo, tetraedro, prisma y pirámide; sin documentar en la guía |
| Desaliasing por sobre-integración | `SPECTRALHPDEALIASING` | Los términos de borde introducen funciones no polinómicas: el aliasing no queda plenamente controlado |
| Disipación de Riemann | `UpwindType` | Ninguno de los diez solvers es entrópicamente estable |

- **La brecha, en una frase verificada** ****: `grep -rn "SVV" solvers/CompressibleFlowSolver/`
  no devuelve nada; los núcleos SVV viven en `library/StdRegions` y son alcanzables sólo por el
  operador de Helmholtz, que el solver compresible nunca ensambla. SVV no está sin usar en el
  compresible: está **inalcanzable**.
- Incluir el diagrama de clases y el punto exacto de inserción.

### 8.3 Punto de inserción

- `ArtificialDiffusion` es una clase base pequeña que ya posee su operador LDG sobre variables
  conservativas y entrega `viscousTensor[j][i] = mu(x) * qfield[j][i]` ****.
  `NonSmoothShockCapture` sólo define el escalar `mu(x)`.
- **SVV entra reemplazando el producto puntual por una operación modal por elemento sobre
  `qfield`**, replicando la proyección que hace `v_SVVLaplacianFilter` pero como filtro autónomo
  del gradiente en lugar de dentro de una operación matricial.
- Consecuencia de alcance: la implementación es una subclase nueva registrada en la factory, con
  cero cambios en el resto del solver.

### 8.4 La arquitectura como objeto de la tesis

*(Decisión tomada: no hay prototipo previo. Se aprende construyendo, y lo construido se explica.)*

- **Por qué no hay prototipo.** El problema técnico que un prototipo debía resolver —el tratamiento
  de trazas— sólo aparece con el término difusivo, es decir con la maquinaria completa; y la parte
  difícil de esa maquinaria ya está resuelta dentro de `ArtificialDiffusion`, que construye su
  operador LDG y entrega el gradiente ****. Un prototipo separado reproduciría trabajo ya hecho y
  agregaría un segundo código que mantener y validar.
- **Qué reemplaza a la función de oráculo del prototipo.** Tres cosas, todas dentro del mismo
  código: los tests unitarios del núcleo y del sensor sobre funciones de prueba con transformada
  modal conocida; el caso de orden de convergencia en flujo suave, que es el detector sensible de
  una inconsistencia en las trazas; y la comparación contra los resultados publicados de los casos
  canónicos (ver `CASOS.md`).
- **Y lo que se gana: la arquitectura pasa de medio a resultado.** El recorrido desde
  `CompressibleFlowSystem::DoOdeRhs` hasta el punto donde se multiplica el gradiente por el
  coeficiente escalar, con el detalle de qué clase hace qué y por dónde viaja cada arreglo, no está
  escrito con ese nivel de detalle en ninguna parte — ni en la guía de usuario ni en la de
  desarrollador. Documentarlo es un aporte con destinatario real, y es lo que hace reproducible la
  implementación.
- Riesgo asumido con esta decisión, para nombrarlo: se aprende C++ y la biblioteca al mismo tiempo,
  sobre el código de producción. La mitigación es que la subclase a escribir es chica y tiene una
  plantilla al lado, `NonSmoothShockCapture` ****.

### 8.5 Implementación en `CompressibleFlowSolver`

- Nueva subclase, registro en la factory, parámetros expuestos en la sesión.
- Cómputo del sensor por elemento; transformación a la base ortogonal; aplicación del núcleo;
  transformación inversa; ensamblado del término difusivo.
- Tratamiento de trazas adoptado (decisión de 7.6).
- Extensión a 2D y 3D: productos tensoriales, elementos deformados, factores geométricos. Nota:
  los núcleos de `StdRegions` cubren todas las formas, incluido el problema de coordenadas
  colapsadas en triángulos **** — ventaja concreta frente al filtro exponencial.
- Costo: dónde se paga (transformadas modales por paso), qué se cachea.

### 8.6 Verificación del código

- Distinguir verificación de *código* (¿resuelve lo que dice resolver?) de validación *física*
  (cap. 9).
- Tests unitarios del núcleo y del sensor.
- Reproducción de los resultados publicados de los casos canónicos — criterio de aceptación. Ver
  `CASOS.md`.
- Orden de convergencia con el operador activo sobre solución suave.
- Conservación de masa, cantidad de movimiento y energía a precisión de máquina. Nota: no afirmar
  que la viscosidad artificial existente es no conservativa — **no hay evidencia de eso** y ambos
  mecanismos se aplican en forma divergencia por operadores DG conservativos ****.
- Diagnósticos nuevos que hay que implementar para el capítulo 10: entropía discreta y su
  derivada, mínimos de densidad y presión, residuo de entropía por elemento. Puntos de partida en
  `library/SolverUtils/Filters/` ****.

### 8.7 Infraestructura de ejecución

- Compilación local y en cluster; dependencias; versión y commit de referencia.
- Salida HDF5, checkpoint/restart, ejecución en paralelo, escalado observado.
- Post-proceso: `FieldConvert`, espectros de energía, integración de fuerzas.
- **Instrumentación de costo desde el día uno**: medición de tiempo de solver excluyendo
  inicialización, post-proceso y E/S, y ejecución de TauBench para normalizar en work units
  ****. Si esto no se instrumenta temprano, el capítulo 9 no tiene números.
- Gestión del conjunto de casos: nomenclatura, trazabilidad de parámetros, scripts de barrido.

---

## Capítulo 9 — `CASOS`

**Función:** definir *antes de correr* qué se corre, con qué malla, qué se mide y qué constituye
éxito.

### 9.1 Estrategia general

- Cinco niveles: verificación de orden → validación canónica → calibración paramétrica → geometría
  canónica con datos experimentales → caso de aplicación.
- Nada sube de nivel sin haber pasado el anterior.

### 9.2 Verificación de orden

- Vórtice isentrópico 2D y advección de onda suave; barrido en `h` y en `P`, con y sin operador.
- Opcionalmente, solución manufacturada de Poiseuille compresible, que es la que usa el propio
  paper del solver (Yan et al., 2021 ****).
- Criterio: el orden observado coincide con el teórico y el operador no lo degrada.

### 9.3 Conjunto mínimo de validación

*(Cada caso está desarrollado por separado en `CASOS.md`: configuración, qué prueba, qué se mide,
referencia y trampas conocidas.)*

**El conjunto de casos se organiza en dos niveles, con propósitos distintos.**

**Nivel 1 — casos que ya trae Nektar++.** Verificado **** contra el árbol de fuentes v5.10.0: el
directorio `solvers/CompressibleFlowSolver/Tests/` contiene 110 casos. Los pertinentes son el vórtice
isentrópico (P1, P3, P8, más variantes con base Gauss y con hexaedro deformado), el flujo de Ringleb
—que trae **solución analítica exacta** implementada en el solver, en régimen transónico—, el
Poiseuille compresible manufacturado —que verifica el **operador difusivo**, que es donde se inserta
el filtro— y el tubo de choque bidimensional en malla mixta, que viene shippeado en **quince
variantes**: viscosidad laplaciana, viscosidad física con sensor de dilatación y con sensor modal,
difusión por penalización interior y LDGNS, y el mismo caso corrido con los ocho solvers de Riemann.

Ventaja: misma malla, mismas condiciones de contorno, mismo flujo de interfaz. **La única variable
que cambia es el mecanismo de estabilización**, que es el experimento controlado que la afirmación de
"igual o mejor" necesita, y el costo de montarlo es casi nulo. Tres resultados salen del tubo de
choque: el contraste directo contra los dos mecanismos existentes; la separación de cuánta disipación
viene del solver de Riemann y cuánta del filtro, usando las ocho variantes ya escritas; y el resultado
de cobertura de forma de elemento, porque esa malla tiene 120 triángulos y el filtro exponencial de
Nektar++ aborta sobre triángulos ****.

**Advertencia que hay que escribir en la tesis:** son **tests de regresión, no de validación**. El
tubo de choque corre un solo paso de `1e-8` con Euler explícito de primer orden y compara normas
contra valores guardados; el vórtice y Ringleb corren diez pasos ****. Son esqueletos con la
geometría y los parámetros ya correctos, que hay que extender a corridas reales.

**Nivel 2 — casos canónicos de la literatura.** Shu–Osher, escalón supersónico y Taylor–Green. Son
necesarios por dos razones: sin un caso compartido con la literatura los resultados no son
contrastables desde afuera del código, y **nada en el test suite de Nektar++ toca turbulencia**
**** — el modo de falla propio del filtro sólo se ve en el Taylor–Green.

*(Respuesta fundada a la pregunta "¿cuál es el mínimo?" — desarrollo en el Anexo A.3.)*

**Núcleo irreducible, cuatro casos.** Es la intersección de las suites de validación de los cinco
trabajos modernos de referencia (Mateo-Gabín et al. 2022; Hennemann et al. 2021; Lin et al. 2023;
Dzanic y Witherden 2022; Yan et al. 2021), todas leídas ****:

1. **Vórtice isentrópico con tabla de orden de convergencia.** Prueba que el estabilizador no
   destruye el orden de diseño donde el flujo es suave. **No negociable.**
2. **Shu–Osher.** Prueba que el sensor se enciende en el choque y se apaga en la región oscilatoria
   suave. Es el caso que discrimina un filtro bien calibrado de uno que borra estructura fina.
3. **Un caso 2D de choque fuerte con estructuras vorticiales**: escalón supersónico M = 3 o doble
   reflexión de Mach. Robustez ante interacción de choques y enrollamiento de capa de corte, sobre
   malla no alineada.
4. **Vórtice de Taylor–Green 3D a Re = 1600.** Prueba que la disipación añadida no corrompe la
   cascada. Es además el caso C3.5 / BS1 de los talleres, con datos de referencia públicos.

**Y, específicamente para una tesis sobre un estabilizador nuevo, tres agregados:**

5. **Un caso de estrés de positividad**: Sod, Sedov o Leblanc.
6. **Demostración del balance de entropía** (ver 7.7), que reemplaza a la demostración analítica
   que esta tesis no puede dar.
7. **Interacción choque–capa límite**, que es lo que usa el paper de referencia del propio solver
   y el caso canónico más cercano a un cohete con capa límite turbulenta.

**Y la condición que hace que el conjunto responda la pregunta de esta tesis:** cada uno de esos
casos debe correrse **tres veces** — DG puro, DG + `NonSmooth`, DG + SVV — a igual costo en work
units. Sin la corrida cabeza a cabeza no hay afirmación de "igual o mejor", sólo de "funciona".

### 9.4 Barrido paramétrico y construcción del mapa

- Variables: amplitud `ε₀`, corte `m_N`, tipo de núcleo (las tres variantes ya implementadas),
  umbral y ancho del sensor, y número de Mach.
- Diseño: no producto cartesiano completo. Barrido en `M` con optimización local de `(ε₀, m_N)`
  alrededor de la predicción analítica de 7.5.
- Métricas por punto: estabilidad (¿termina la corrida?), error contra referencia, disipación del
  operador frente a disipación física, y costo.
- Producto: la función `f(M)` y su banda de admisibilidad.

### 9.5 Caso tridimensional canónico: Taylor–Green

- Setup estándar, de la especificación del taller (Hillewaert, caso C3.5) ****: dominio
  periódico `[−πL, πL]³`, `Re = 1600`, `M₀ = 0.10`, `γ = 1.4`, `Pr = 0.71`, viscosidad de volumen
  nula, hasta `t = 20 t_c` con pico de disipación en `t ≈ 8 t_c`, malla base ≈ 256³ grados de
  libertad.
- Métricas obligatorias del taller ****: energía cinética, tasa de disipación de energía
  cinética, enstrofía, y norma de vorticidad sobre la cara `x/L = −π` en `t/t_c = 8`.
- **Métrica adicional que esta tesis necesita** ****: la descomposición de la disipación en
  `ε₁` (deviatórica), `ε₂` (viscosidad de volumen, nula acá) y `ε₃ = −(1/ρ₀Ω)∫ p ∇·v dΩ`, que es
  la dilatación de presión. Reportar `ε₁` y `ε₃` por separado es la forma sancionada por el taller
  de mostrar cómo entra la compresibilidad, y conecta este capítulo directamente con el 5.
- El taller también sugiere reportar la variante numérica de la tasa de disipación, incluyendo los
  términos de salto del esquema DG, y compararla con la consistente ****. **Ése es exactamente
  el diagnóstico que cuantifica cuánta disipación agrega el operador.**
- Datos de referencia: solución pseudo-espectral del taller convergida en 512³ ****;
  van Rees et al. (2011), *JCP* 230:2794–2805 ****; Brachet et al. (1983), *JFM* 130:411–452
  **** (seis autores, no cinco; el rango de páginas es 411–452).
- **Caso puente hacia el extremo supersónico:** el TGV compresible a M = 1.25 de Lusher y Sandham
  (2021), *AIAA J.* 59:533–545 **[V, sólo resumen]**, es el único TGV donde la disipación
  dilatacional y los *shocklets* son relevantes. Es el puente natural entre la validación a bajo
  Mach y el cohete, y es el régimen donde SVV debe hacer las dos tareas a la vez. **Verificar la
  condición inicial y los datos de referencia en el paper antes de adoptarlo.**

### 9.6 Escalón intermedio: HB-2

*(Recomendación nueva en v0.3.)*

- **Qué es:** modelo estándar hipersónico HB-2, cono romo de semiángulo 25° + cilindro + flare de
  10°, longitud 4.9 D. Topológicamente, un cohete sin aletas.
- **Por qué conviene:** hay datos experimentales de fuerza publicados de M = 1.5 a 10 (Gray y
  Lindsay, AEDC-TDR-63-137; Gray, AEDC-TDR-64-137 ****), extendidos al rango transónico
  M = 0.7–1.4 por Damljanović et al. (2025), *Aerospace* 12:131 ****. Y —esto es lo decisivo—
  las tablas separan **coeficiente axial de antecuerpo `C_Af` y coeficiente de presión de base
  `C_pb`**, que es exactamente la descomposición que hace falta para el cohete.
- **Función dentro del conjunto:** si la comparación con el Aconcagua da 17 % de discrepancia, HB-2 es
  lo que dice si el problema es el solver o la telemetría. Sin ese escalón, la discrepancia no es
  interpretable.
- Es además una geometría axisimétrica que una malla curvada DG maneja limpiamente.

### 9.7 Caso de aplicación: cohete Aconcagua

- **Descripción.** Vehículo supersónico desarrollado en la facultad; geometría en STL; telemetría
  de vuelo disponible. Justificación de la elección: geometría real, régimen dentro del rango, y
  una comunidad concreta que usa el resultado.
- **Preparación geométrica.** Limpieza del STL, simplificaciones declaradas (aletas, protuberancias,
  tobera), dominio de cálculo y condiciones de contorno.
- **Malla de alto orden.** NekMesh (Green et al., 2024, *CPC* 298:109089 ****): malla lineal
  de base, curvado de superficie contra el STL, capa límite, y **verificación de jacobiano positivo
  en los elementos curvados** — el punto de falla más frecuente de esta etapa, y reconocido como
  tal en la literatura: Wang et al. (2013) **** listan la generación de mallas de alto orden como
  el primer *pacing item* y señalan que *"the main difficulty is that cells near the curved
  geometries can overlap each other."*
- **Malla de OpenFOAM.** Generada desde el mismo STL. Precedente citable para el desajuste: Coppeans
  et al. (2023) **** resuelven el problema degradando la malla curvada a lineal conectando los
  nodos de superficie de alto orden. Documentar que el código FVM ve una aproximación facetada de
  la misma superficie.
- **Condiciones de vuelo.** Puntos de la trayectoria a lo largo de **todo el vuelo**, declarando
  `M`, altitud, `Re` y ángulo de ataque supuesto de cada uno, y señalando cuáles caen en tramo
  propulsado y cuáles en coasteo — porque la incertidumbre de la referencia no es la misma en los
  dos (3.6). En los puntos propulsados, declarar la curva de empuje supuesta.
- **Cantidades de interés.** Coeficiente de arrastre (métrica primaria, porque es lo que la
  telemetría entrega), distribución de presión, posición del choque de proa, presión de base.
- **Protocolo de comparación.** Ver 9.8.
- **Presupuesto de incertidumbre de la referencia experimental.** Enumerar y, donde se pueda,
  cuantificar las fuentes de 3.7. Esto no es relleno: es lo que separa "el CFD dio 8 % distinto"
  de "el CFD está dentro del margen de error del dato".

### 9.8 Protocolo de comparación de costo

*(Respuesta fundada a "¿no deberíamos partir de una malla igual?" — desarrollo en el Anexo A.4.)*

- **Base principal: error frente a costo.** Wang et al. (2013) **** rechazan explícitamente
  tanto la malla igual como los grados de libertad iguales como base de una afirmación de
  eficiencia. Se reportan curvas error–work units para cada solver, variando resolución (`h` y `P`
  en Nektar++, `h` en OpenFOAM).
- **Base secundaria: error frente a grados de libertad.** Independiente del hardware. Los talleres
  definen `h = N_DOF^(−1/d)` justamente para que esa curva sea comparable entre métodos ****.
- **Malla igual: sí, pero como experimento de control etiquetado, en apéndice.** La aritmética que
  lo obliga: un hexaedro `P = 4` con base tensorial lleva `(P+1)³ = 125` puntos de solución; una
  celda FV de segundo orden lleva 1. "La misma malla de 100 mil celdas" le da al solver DG 12.5
  millones de grados de libertad contra 100 mil — una ventaja de resolución de 125× disfrazada de
  equidad. Vermeire et al. (2017) **** lo enuncian operativamente: *"we use fewer elements in
  the PyFR mesh since the FR scheme has multiple solution points per element."*
  Lo que la malla igual sí sirve: como prueba de robustez ("¿el solver DG sobrevive sobre la malla
  del FVM?") y como terreno común honesto cuando la geometría es fija. Reportarlo así, con el
  descargo citado.
- **Definición de costo** ****: work units = `NP · T_solver / T_TauBench`, excluyendo
  inicialización, post-proceso y E/S. El tiempo de generación de malla **no** entra —ninguna
  comparación publicada lo incluye— pero sí se discute cualitativamente, y esa discusión honesta
  suma en lugar de restar.
- **Regla del 30 %** ****: no afirmar diferencia de eficiencia menor a ese umbral.
- **Igualar condiciones de ejecución**: misma máquina, mismo número de procesos, mismo criterio de
  convergencia o mismo tiempo físico simulado, misma tolerancia.
- **Amenazas específicas**: madurez desigual de implementaciones, esfuerzo de ajuste desigual,
  escalado en paralelo distinto. Precedente y advertencia: la comparación Nektar++ contra OpenFOAM
  de Jiang y Cheng (2021) reporta un orden de magnitud a favor de Nektar++, pero la cifra de
  OpenFOAM es **extrapolada** ****. No repetir ese atajo.

### 9.9 Reproducibilidad

- Versiones, commits, archivos de sesión, semillas, scripts de post-proceso, ubicación de los datos.

---

# PARTE 4 — CIERRE

## Capítulo 10 — `RESULTADOS`

**Regla:** ordenar por afirmación, no por caso. Cada sección responde una sub-pregunta de 3.1.

### 10.1 Verificación del operador

- Órdenes de convergencia con y sin operador; conservación; comportamiento del sensor; resultado
  de los tests unitarios del núcleo.

### 10.2 Estabilización ante discontinuidades

- Resultados 1D y 2D: perfiles, oscilación residual, ancho del choque en número de elementos.
- **Comparación cabeza a cabeza contra `NonSmooth` a igual costo.** Es la sección que responde
  "¿igual o mejor?".
- Comparación entre las tres variantes de núcleo.

### 10.3 Régimen turbulento sub-resuelto

- Taylor–Green: energía cinética, tasa de disipación, enstrofía, espectros contra `k^(−5/3)` con
  el número de onda de Nyquist marcado.
- Descomposición `ε₁` / `ε₃` y variante numérica de la disipación con términos de salto.
- Verificación del balance disipación numérica frente a física.

### 10.4 Balance de entropía y admisibilidad

- `dS/dt` para DG puro, DG + AV existente y DG + SVV.
- Mínimos de densidad y presión; catálogo de violaciones de positividad observadas.

### 10.5 Mapa paramétrico frente al número de Mach

- La función `f(M)` y su banda de admisibilidad.
- Contraste contra la predicción analítica de 7.5: dónde acierta y dónde no.
- Identificación de transiciones de régimen, si las hay.

### 10.6 HB-2

- Coeficiente axial de antecuerpo y coeficiente de presión de base contra datos experimentales, en
  el rango M = 1.5–3.
- Es el resultado que calibra la interpretación de 10.7.

### 10.7 Aconcagua

- Campos y cantidades de interés en cada punto de vuelo, separando tramo propulsado y coasteo.
- `C_D(M)` de CFD contra `C_D(M)` de telemetría, con barras de incertidumbre de ambos lados.
- Curvas error–costo de ambos solvers; costo a error fijo.

### 10.8 Fallos y casos no resueltos

- Configuraciones que divergieron, fallos de positividad, límites del sensor, casos recortados.
- **Sección obligatoria:** su ausencia es lo primero que un jurado nota.

---

## Capítulo 11 — `DISCUSION`

- **Interpretación del mapa paramétrico.** ¿Qué explica la forma de `f(M)`? Conectar con el
  capítulo 5. Hipótesis a contrastar: si Duan et al. tienen razón y los términos dilatacionales son
  pequeños en capa adherida, entonces `f(M)` debería estar dominada por la escala de velocidad
  característica local `|u| + a` y no por efectos de compresibilidad turbulenta. Si los datos dicen
  otra cosa, es numérico, no físico — y hay que investigarlo.
- **¿Se cumplieron los criterios de trabajo de 3.4?** Uno por uno, con el resultado a la vista.
- **Lectura de la comparación de costo.** Qué dice y qué no: el resultado es sobre *esta*
  implementación, *esta* geometría y *este* rango. Cuantificar cuánto del costo es intrínseco al
  método y cuánto es falta de optimización, usando el perfilado de 3.7.
- **Comparación con la literatura.** Contraste con Mateo-Gabín et al. (2022), Manzanero et al.
  (2020) y Kirby y Sherwin (2006). Punto central: esos trabajos corren a Mach ≈ 0.1 y sobre
  esquemas de forma split; esta tesis corre hasta M = 3 sobre DG estándar. Qué reproduce, qué
  extiende, dónde difiere y por qué.
- **Limitaciones.** Ausencia de garantía de estabilidad entrópica del esquema completo, ausencia de
  limitador de positividad, rango de `P` explorado, dependencia del sensor de la variable elegida,
  y el hecho de que la validación del cohete no es una LES validada.
- **Implicaciones prácticas.** Recomendación operativa: qué valores usar como punto de partida para
  un caso nuevo en el rango, y qué verificar antes de confiar en el resultado.

---

## Capítulo 12 — `CONCLUSION`

- Respuesta directa a la pregunta de 3.1, en un párrafo y sin hedging.
- Aportes logrados, enfrentados uno a uno con la lista de 1.4, incluyendo los parciales.
- Contribución al problema de acceso de 1.1: qué puede hacer hoy un equipo pequeño que antes no
  podía, y a qué costo.
- Trabajo futuro, ordenado por madurez:
  1. Viscosidad artificial de corrección de entropía sobre DG estándar, siguiendo Chan (2025) —
     arquitectónicamente la misma forma que el operador ya implementado.
  2. Limitador de positividad tipo Zhang–Shu.
  3. Discretización de volumen en forma split con SBP, que abriría la ruta de Mateo-Gabín y de
     Lin, Chan y Tomas. Es un proyecto en sí mismo.
  4. Extensión a M > 3 y a efectos de gas real.
  5. Optimización de la implementación y escalado en cluster.
  6. Adaptatividad `p` guiada por el mismo sensor modal — la infraestructura ya existe en el
     `DriverAdaptive` de Nektar++ ****.
- Cierre: volver al Aconcagua. La tesis empezó con un equipo que necesita coeficientes confiables
  sin cómputo masivo; decir si ahora los tiene.

---
# Bibliografía verificada

*Ordenada por rol en la tesis. El estado de verificación es el de agosto de 2026 y corresponde al
trabajo de relevamiento hecho para este documento, no a lectura completa por el autor. Todo lo
marcado **[F]** o **[X]** debe leerse antes de citarse.*

**Fundamento del operador**

1. Tadmor, E. (1989). Convergence of spectral methods for nonlinear conservation laws.
   *SIAM J. Numer. Anal.* 26(1):30–44. **** — PDF local.
2. Maday, Y., Ould Kaber, S.M., Tadmor, E. (1993). Legendre pseudospectral viscosity method for
   nonlinear conservation laws. *SIAM J. Numer. Anal.* 30(2):321–342. **** — PDF local.
3. Kirby, R.M., Sherwin, S.J. (2006). Stabilisation of spectral/hp element methods through spectral
   vanishing viscosity. *CMAME* 195:3128–3144, doi:10.1016/j.cma.2004.09.019. **** — PDF local.
   **Incompresible.**
4. Moura, R.C., Sherwin, S.J., Peiró, J. (2016). Eigensolution analysis of spectral/hp continuous
   Galerkin approximations to advection–diffusion problems: insights into spectral vanishing
   viscosity. *JCP* 307:401–422, doi:10.1016/j.jcp.2015.12.009. **** (resumen).
5. Karamanos, G.-S., Karniadakis, G.E. (2000). A spectral vanishing viscosity method for large-eddy
   simulations. *JCP* 163(1):22–50, doi:10.1006/jcph.2000.6552. ****, texto completo **[X]**.
   **Leer antes de citar.** *Nota: Kirby y Sherwin (2006) lo citan como volumen 162; el correcto es
   163.*
6. Pasquetti, R. (2006). Spectral vanishing viscosity method for LES: sensitivity to the SVV
   control parameters. *J. Sci. Comput.* 27(1–3):365–375, doi:10.1007/s10915-005-9029-9. **** —
   PDF local. **Incompresible.**

**Sensor y captura de choques**

7. Persson, P.-O., Peraire, J. (2006). Sub-cell shock capturing for discontinuous Galerkin methods.
   AIAA 2006-112, doi:10.2514/6.2006-112. **** — PDF local.
8. Hennemann, S., Rueda-Ramírez, A.M., Hindenlang, F.J., Gassner, G.J. (2021). A provably entropy
   stable subcell shock capturing approach for high order split form DG. *JCP* 426:109935. ****.

**SVV en DG compresible**

9. Manzanero, J., Ferrer, E., Rubio, G., Valero, E. (2020). Design of a Smagorinsky spectral
   vanishing viscosity turbulence model for discontinuous Galerkin methods. *Computers & Fluids*
   200:104440, doi:10.1016/j.compfluid.2020.104440. **** — manuscrito abierto en oa.upm.es.
   **Corre a M = 0.1 en todos los casos.**
10. Mateo-Gabín, A., Manzanero, J., Valero, E. (2022). An entropy stable spectral vanishing
    viscosity for discontinuous Galerkin schemes: application to shock capturing and LES models.
    *JCP* 471:111618, doi:10.1016/j.jcp.2022.111618. **** — PDF local (archivo mal nombrado).
    **Tres autores. Corregir la entrada previa que listaba a Rueda-Ramírez y Rubio.**

**Estabilidad entrópica y positividad**

11. Lin, Y., Chan, J., Tomas, I. (2023). A positivity preserving strategy for entropy stable
    discontinuous Galerkin discretizations of the compressible Euler and Navier-Stokes equations.
    *JCP* 475:111850, doi:10.1016/j.jcp.2022.111850. ****.
12. Chan, J. (2025). An artificial viscosity approach to high order entropy stable discontinuous
    Galerkin methods. *JCP* 543:114380, doi:10.1016/j.jcp.2025.114380. ****. **La ruta futura.**
13. Chan, J. (2018). On discretely entropy conservative and entropy stable discontinuous Galerkin
    methods. *JCP* 362:346–374, doi:10.1016/j.jcp.2018.02.033. ****.
14. Gassner, G.J., Winters, A.R., Kopriva, D.A. (2016). Split form nodal discontinuous Galerkin
    schemes with summation-by-parts property for the compressible Euler equations. *JCP*
    327:39–66. ****.
15. Zhang, X., Shu, C.-W. (2010). On maximum-principle-satisfying high order schemes for scalar
    conservation laws / positivity-preserving schemes for compressible Euler. *JCP* 229(23):8918–8934.
    ****. Y (2011), *Proc. R. Soc. A* 467:2752–2776. ****.

**Turbulencia y escalas**

16. Beck, A.D., Bolemann, T., Flad, D., Frank, H., Gassner, G.J., Hindenlang, F., Munz, C.-D.
    (2014). High-order discontinuous Galerkin spectral element methods for transitional and
    turbulent flow simulations. *IJNMF* 76(8):522–548, doi:10.1002/fld.3943. **** (resumen).
17. Moura, R.C., Sherwin, S.J., Peiró, J. (2015). Linear dispersion–diffusion analysis and its
    application to under-resolved turbulence simulations using DG spectral/hp methods. *JCP*
    298:695–710, doi:10.1016/j.jcp.2015.06.020. **** (resumen).
18. Moura, R.C., Mengaldo, G., Peiró, J., Sherwin, S.J. (2017). On the eddy-resolving capability of
    high-order DG approaches to implicit LES / under-resolved DNS of Euler turbulence. *JCP*
    330:615–623, doi:10.1016/j.jcp.2016.10.056. **** (resumen). *La "regla del 1 %".*
19. Fernandez, P., Moura, R.C., Mengaldo, G., Peraire, J. (2019). Non-modal analysis of spectral
    element methods: towards accurate and robust large-eddy simulations. *CMAME* 346:43–62.
    **** (resumen).
20. Rubio, G., Ntoukas, G., Chávez-Módena, M., Mariño, O.A., Font, B., Lehmkuhl, O., Valero, E.,
    Ferrer, E. (2025). Can explicit subgrid models enhance implicit LES simulations? arXiv:2512.04574.
    ****. *Fuente de la divergencia de DG estándar sobre TGV.*
21. Winters, A.R., Moura, R.C., Mengaldo, G., Gassner, G.J., Walch, S., Peiró, J., Sherwin, S.J.
    (2018). A comparative study on polynomial dealiasing and split form DG schemes for
    under-resolved turbulence computations. *JCP* 372:1–21. **** (resumen).
22. Duan, L., Beekman, I., Martín, M.P. (2011). DNS of hypersonic turbulent boundary layers.
    Part 3. Effect of Mach number. *JFM* 672:245–267, doi:10.1017/S0022112010005902. ****
    (resumen). *La justificación de despreciar términos dilatacionales en capa adherida.*
23. Sarkar, S., Erlebacher, G., Hussaini, M.Y., Kreiss, H.O. (1991). The analysis and modelling of
    dilatational terms in compressible turbulence. *JFM* 227:473–493. **** — PDF local.
24. Lele, S.K. (1994). Compressibility effects on turbulence. *Annu. Rev. Fluid Mech.* 26:211–254.
    ****.
25. Erlebacher, G., Hussaini, M.Y., Speziale, C.G., Zang, T.A. (1992). Toward the LES of
    compressible turbulent flows. *JFM* 238:155–185. ****.
26. Garnier, E., Adams, N., Sagaut, P. (2009). *Large Eddy Simulation for Compressible Flows*.
    Springer, doi:10.1007/978-90-481-2819-8. ****, contenido **[F]**.

**Nektar++**

27. Cantwell, C.D. et al. (2015). Nektar++: an open-source spectral/hp element framework. *CPC*
    192:205–219. **** — PDF local.
28. Moxey, D. et al. (2020). Nektar++: enhancing the capability and application of high-fidelity
    spectral/hp element methods. *CPC* 249:107110 (arXiv:1906.03489). **** (preprint).
    *Verificar números de página de la versión publicada.*
29. Yan, Z.-G., Pan, Y., Castiglioni, G., Hillewaert, K., Peiró, J., Moxey, D., Sherwin, S.J.
    (2021). Nektar++: design and implementation of an implicit, spectral/hp element, compressible
    flow solver using a Jacobian-free Newton Krylov approach. *Comput. Math. Appl.* 81:351–372
    (arXiv:2002.04222). **** (preprint). *Suite de validación de referencia.*
30. Green, M.D. et al. (2024). NekMesh: an open-source high-order mesh generation framework.
    *CPC* 298:109089. ****.

**Metodología de comparación**

31. Wang, Z.J. et al. (2013). High-order CFD methods: current status and perspective. *IJNMF*
    72(8):811–845, doi:10.1002/fld.3767. ****. *La fuente de la metodología error–costo y de la
    regla del 30 %.*
32. Guías de los talleres HiOCFD4 (2016) y HiOCFD5 (2017), how4.cenaero.be y how5.cenaero.be.
    ****. *Definición de work units y de `h = N_DOF^(−1/d)`.*
33. Hillewaert, K. (2011). Problem C3.5: DNS of the Taylor–Green vortex at Re = 1600.
    1st International Workshop on High-Order CFD Methods. ****. *Especificación y datos de
    referencia a 512³.*
34. Vermeire, B.C., Witherden, F.D., Vincent, P.E. (2017). On the utility of GPU accelerated
    high-order methods for unsteady flow simulations: a comparison with industry-standard tools.
    *JCP* 334:497–521. ****.
35. Coppeans, A.W., Fidkowski, K.J., Martins, J.R.R.A. (2023). Comparison of finite volume and
    high-order DG based aerodynamic shape optimization. AIAA 2023-1845. **** (parcial).
36. Jiang, H., Cheng, L. (2021). Large-eddy simulation of flow past a circular cylinder for
    Reynolds numbers 400 to 3900. *Phys. Fluids* 33(3), doi:10.1063/5.0041168. ****;
    cifras de costo desde el sitio de Nektar++ ****. *La cifra de OpenFOAM es extrapolada.*

**Casos de validación y referencia experimental**

37. Brachet, M.E., Meiron, D.I., Orszag, S.A., Nickel, B.G., Morf, R.H., Frisch, U. (1983).
    Small-scale structure of the Taylor–Green vortex. *JFM* 130:411–452. ****. *Seis autores;
    páginas 411–452.*
38. van Rees, W.M., Leonard, A., Pullin, D.I., Koumoutsakos, P. (2011). *JCP* 230(7):2794–2805.
    ****.
39. Lusher, D.J., Sandham, N.D. (2021). Assessment of low-dissipative shock-capturing schemes for
    the compressible Taylor–Green vortex. *AIAA J.* 59(2):533–545, doi:10.2514/1.J059672. ****
    (resumen). *El TGV a M = 1.25. Verificar condición inicial antes de adoptarlo.*
40. Damljanović, D., Vuković, Đ., Ocokoljić, G., Rašuo, B. (2025). New transonic tests of HB-2
    hypersonic standard models in the VTI T-38 trisonic wind tunnel. *Aerospace* 12(2):131,
    doi:10.3390/aerospace12020131. ****, acceso abierto.
41. Gray, J.D., Lindsay, E.E. (1963). Force tests of standard hypervelocity ballistic models HB-1
    and HB-2 at Mach 1.5 to 10. AEDC-TDR-63-137. ****. Y Gray, J.D. (1964), AEDC-TDR-64-137.
    ****.
42. Bawa, V., Michalski, Q. (2025). An experimental framework for in-flight determination of rocket
    aerodynamic drag. AIAA 2025-104785, doi:10.2514/6.2025-104785. ****.
43. Kierulf, B. (2022). Critical review of sounding rocket flight data and computational fluid
    dynamics efforts. ICAS 2022, paper 0414. ****.
44. Slotnick, J. et al. (2014). CFD Vision 2030 study: a path to revolutionary computational
    aerosciences. NASA/CR-2014-218178. ****.
45. Statnikov, V., Sayadi, T., Meinke, M., Schmid, P., Schröder, W. (2015). *Phys. Fluids*
    27:016103. ****. Y Simon, F., Deck, S., Guillen, P., Sagaut, P. (2006), *AIAA J.*
    44(11):2578–2590. ****. *Práctica de resolución de escalas en lanzadores.*

**Marcadas para verificar antes de usar**

- Gassner, G.J., Beck, A.D. (2013). *TCFD* 27(3–4):221–237. DOI verificado; sin resumen accesible;
  texto no leído. **No citar afirmaciones puntuales.**
- Karniadakis, G.E., Sherwin, S.J. (2005). *Spectral/hp Element Methods for CFD*, 2ª ed., OUP.
  Sin acceso; ninguna cita textual disponible.
- Chapelier, J.-B., de la Llave Plata, M., Lamballais, E. (2016). *CMAME* 307:275–299. **[F]**.
- Vreman, B., Geurts, B., Kuerten, H. (1995). *Appl. Sci. Res.* 54:191–203. ****, no leído.
- AGARD-AR-303 y AGARD-AR-138: registros confirmados, listas de casos **no** verificadas.
- Morkovin (1962): capítulo de libro sin DOI, rango de páginas sin confirmar.
- No existe, hasta donde se pudo verificar, ningún trabajo de Chapelier sobre SVV en DG compresible.
  **Tratar como inexistente hasta encontrarlo.**
