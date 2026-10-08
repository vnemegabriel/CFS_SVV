# Integración en el build y estrategia de testeo

Referencias de línea verificadas contra tu copia en `/opt/nektar++`.

---

## 1. Dónde van los archivos

```
solvers/CompressibleFlowSolver/
├── SVV/
│   ├── SVVOperator.h        SVVOperator.cpp
│   └── SVVExpKernel.h       SVVExpKernel.cpp
└── Turbulence/
    ├── TurbulenceModel.h    TurbulenceModel.cpp
    └── Smagorinsky.h        Smagorinsky.cpp
```

Sigue la convención de las carpetas hermanas (`ArtificialDiffusion/`,
`RiemannSolvers/`, `BoundaryConditions/`): una carpeta por familia de clases,
con la base y las derivadas juntas.

---

## 2. Cambios en `solvers/CompressibleFlowSolver/CMakeLists.txt`

Dos ediciones, ambas alfabéticas.

**(a) Fuentes** — dentro de `SET(CompressibleFlowSolverSource ...)`, que hoy
arranca en la línea 10. Agregar antes del cierre `)` de la línea 43:

```cmake
        SVV/SVVOperator.cpp
        SVV/SVVExpKernel.cpp
        Turbulence/TurbulenceModel.cpp
        Turbulence/Smagorinsky.cpp
```

**(b) Headers** — dentro de `SET(CompressibleFlowSolverHeaders ...)`, cerca de
la línea 71 donde están los de `ArtificialDiffusion`:

```cmake
        SVV/SVVOperator.h
        SVV/SVVExpKernel.h
        Turbulence/TurbulenceModel.h
        Turbulence/Smagorinsky.h
```

Después de editar, `cmake` se re-ejecuta solo en el próximo `make`. Si no lo
hace, borrá `CMakeCache.txt` — no borres todo el directorio de build.

**Trampa a anticipar:** si el ejecutable compila pero tu clase "no existe" en
tiempo de ejecución (la factory dice que no conoce la clave), es casi seguro
que el linker descartó el objeto porque nada referencia su símbolo. Se
resuelve asegurando que el `.cpp` esté en `SOURCES` y no solo en
`LIBRARY_SOURCES`.

---

## 3. Puntos de enganche exactos

### 3.1 SVV

Archivo: `EquationSystems/NavierStokesCFE.cpp` (y su gemelo implícito
`NavierStokesImplicitCFE.cpp` si vas a soportar los dos).

1. **Miembro** en `NavierStokesCFE.h`, junto a los existentes:
   ```cpp
   SVVOperatorSharedPtr m_svvOperator;
   ```

2. **Construcción** en `NavierStokesCFE::v_InitObject()`, imitando cómo se
   construye `m_artificialDiffusion` en `CompressibleFlowSystem`:
   ```cpp
   if (m_session->DefinesSolverInfo("SVVType"))
   {
       std::string svvType = m_session->GetSolverInfo("SVVType");
       m_svvOperator = GetSVVOperatorFactory().CreateInstance(
           svvType, m_session, m_fields, m_spacedim);
   }
   ```

3. **Aplicación** en `NavierStokesCFE::v_DoDiffusion` (declarado en
   `NavierStokesCFE.h:128`), después del término viscoso físico:
   ```cpp
   if (m_svvOperator)
   {
       m_svvOperator->DoSVV(inarray, outarray);
   }
   ```
   El `if (m_svvOperator)` funciona porque un `shared_ptr` vacío convierte a
   `false`. Es el idioma estándar del codebase para operadores opcionales.

### 3.2 Turbulencia

El enganche es más quirúrgico y por eso más seguro.

Archivo: `EquationSystems/NavierStokesCFE.cpp`, función
`GetViscosityAndThermalCondFromTemp` (declarada en `NavierStokesCFE.h:155`).
Hoy hace, en esencia:

```cpp
mu = m_varConv->GetDynamicViscosity(temperature);   // NavierStokesCFE.h:192
thermalCond = cp / Pr * mu;
```

Pasa a:

```cpp
mu = m_varConv->GetDynamicViscosity(temperature);
thermalCond = cp / Pr * mu;

if (m_turbulenceModel)
{
    Array<OneD, NekDouble> muT(nPts, 0.0);
    Array<OneD, NekDouble> kappaT(nPts, 0.0);

    m_turbulenceModel->GetEddyViscosity(physfield, derivatives, muT);
    m_turbulenceModel->GetTurbulentThermalConductivity(muT, kappaT);

    Vmath::Vadd(nPts, mu, 1, muT, 1, mu, 1);
    Vmath::Vadd(nPts, thermalCond, 1, kappaT, 1, thermalCond, 1);
}
```

**Complicación real a resolver:** `GetViscosityAndThermalCondFromTemp` recibe
solo la temperatura, no los gradientes de velocidad. Vas a tener que:

- o bien mover el cálculo de `mu_t` a `GetViscousFluxVector`, donde los
  gradientes sí están disponibles;
- o bien guardar `mu_t` como miembro, calculado una vez por paso justo antes
  (patrón que ya usa el solver para la viscosidad artificial).

La segunda opción es menos invasiva. Decidilo mirando el orden real de
llamadas — es exactamente el ejercicio de trazado de la Fase 2 del plan.

---

## 4. XML de sesión

```xml
<CONDITIONS>
  <SOLVERINFO>
    <I PROPERTY="EQType"          VALUE="NavierStokesCFE" />
    <I PROPERTY="SVVType"         VALUE="ExpKernel"       />
    <I PROPERTY="SVVMomentumOnly" VALUE="True"            />
    <I PROPERTY="TurbulenceModel" VALUE="Smagorinsky"     />
    <I PROPERTY="SGSShockLimiter" VALUE="True"            />
  </SOLVERINFO>

  <PARAMETERS>
    <P> SVVCutoffRatio = 0.75 </P>
    <P> SVVDiffCoeff   = 0.1  </P>
    <P> Cs             = 0.1  </P>
    <P> PrT            = 0.9  </P>
  </PARAMETERS>
</CONDITIONS>
```

Los nombres de parámetros son deliberadamente los mismos que usa el solver
incompresible (`SVVCutoffRatio`, `SVVDiffCoeff` — ver
`VelocityCorrectionScheme.cpp:1110` y `:1064`). Reusar nomenclatura le ahorra
trabajo a quien venga después y le simplifica la vida a vos si algún día
querés mandar un merge request upstream.

---

## 5. Escalera de tests

### Nivel 1 — Unit tests del kernel

Ubicación: `library/UnitTests/` (mirá los existentes para el formato Boost.Test).

Para `SVVExpKernel::GetKernel`, tres propiedades que no requieren correr CFD:

| Propiedad | Chequeo |
|---|---|
| Los modos resueltos no se tocan | `kernel[n] == 0.0` exacto para todo `n <= M` |
| El modo más alto recibe disipación plena | `kernel[P] == 1.0` |
| Transición monótona | `kernel[n+1] > kernel[n]` para `n > M` |

Para `GetDeviatoricStrainRate`, campos analíticos:

| Campo | `|S|` esperado |
|---|---|
| Corte puro `u = (y, 0, 0)` | 1.0 |
| Compresión pura `u = (x, y, z)` | **0.0** ← este es el test que atrapa el bug de no restar la traza |
| Rotación rígida `u = (-y, x, 0)` | 0.0 |

El segundo es el más importante de los tres: si te da distinto de cero, tu
modelo va a meter viscosidad turbulenta en los choques.

### Nivel 2 — Regresión

Copiá un `.tst` existente. El formato es (ver
`Tests/CylinderSubsonic_P3.tst`):

```xml
<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Navier-Stokes, cilindro subsonico, SVV ExpKernel P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_SVV_P3.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_SVV_P3.xml</file>
        <file description="Restart File">CylinderSubsonic_SVV_P3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-12">...</value>
            <value variable="rhou" tolerance="1e-12">...</value>
            <value variable="rhov" tolerance="1e-12">...</value>
            <value variable="E"    tolerance="1e-12">...</value>
        </metric>
    </metrics>
</test>
```

Los valores se obtienen corriendo el caso una vez, verificando **a mano** que
el resultado es físicamente sensato, y recién ahí congelándolos. Un test de
regresión no valida que el resultado sea correcto — valida que no cambió. No
lo confundas con verificación.

Correr: `ctest -R CylinderSubsonic_SVV` desde el directorio de build.

### Nivel 3 — Consistencia (el test más barato y más útil)

Con `SVVDiffCoeff = 0` y `Cs = 0`, el solver debe reproducir el resultado
original **bit a bit**. Corré el mismo caso con y sin tus clases activas y
diferenciá los `.fld`. Si difieren, tenés un bug de plomería (índice mal,
array sin inicializar, suma en vez de sobreescritura) y lo encontrás en
minutos en vez de en semanas.

Hacé este test después de **cada** cambio estructural.

### Nivel 4 — Orden de convergencia

Vórtice isentrópico (ya está implementado: `BoundaryConditions/IsentropicVortexBC`).
Barré P = 2..8 con el mismo caso, graficá error L2 vs P en escala semilog.

- Sin SVV: la pendiente debe ser exponencial (convergencia espectral).
- Con SVV y coeficiente chico: la pendiente debe degradarse **poco**. Si se
  degrada mucho, tu cutoff está demasiado bajo y estás disipando modos
  resueltos — que es exactamente lo que SVV promete no hacer.

Este gráfico va sí o sí en el escrito: es la justificación cuantitativa de por
qué SVV y no viscosidad artificial.

### Nivel 5 — Física, subsónico

Taylor–Green a Re=1600. **Ya tenés el caso montado** en
`code/3d_mortensen_dns/` con `tgv_final.h5` y `energy_spectrum.png`.

Dos métricas:

1. **Tasa de disipación de energía cinética** `-dk/dt` vs. tiempo, contra los
   datos DNS de referencia de van Rees et al. / Wang et al. El pico alrededor
   de t≈9 es la firma del caso.
2. **Espectro de energía**: verificar la pendiente −5/3 en el rango inercial y
   observar qué le hace SVV al extremo de alto número de onda. El efecto debe
   ser una caída limpia solo en los últimos modos, sin tocar el rango inercial.

Barrido mínimo a reportar: `SVVDiffCoeff` ∈ {0, 0.1, 0.5, 1.0} y
`SVVCutoffRatio` ∈ {0.5, 0.75, 0.9}.

### Nivel 6 — Física, supersónico

1. **Turbulencia isótropa compresible decayente** (configuración tipo Samtaney
   et al.), barriendo Mach turbulento M_t ∈ {0.1, 0.3, 0.6}. Es el caso canónico
   para verificar que un modelo de submalla compresible se comporta al variar
   la compresibilidad.
2. **Un caso con choque limpio** (tubo de choque, o cuña oblicua) sin
   turbulencia. Objetivo: verificar que el modelo LES **no** altera el salto de
   Rankine–Hugoniot ni ensancha el choque. Si lo hace, el limitador de choque
   no está funcionando.
3. **Descomposición de disipación.** Graficá en el mismo campo:
   `mu_laminar`, `mu_artificial` (captura de choque), `mu_SVV` y `mu_t`.
   Es para esto que las clases exponen sus viscosidades por separado. Este
   gráfico es probablemente el resultado más defendible de toda la tesis: es
   la respuesta directa a la pregunta "¿cómo sabés que tu modelo está haciendo
   algo y no es solo disipación numérica?".

---

## 6. Orden de trabajo recomendado

| Paso | Qué | Test que lo cierra |
|---|---|---|
| 1 | Compilar Nektar++ y correr un caso existente | Corre sin error |
| 2 | Agregar las clases al build con los cuerpos vacíos | Compila y linkea |
| 3 | Verificar que la factory las encuentra desde el XML | Un `cout` en el constructor se imprime |
| 4 | Implementar `v_GetKernel` | Nivel 1 |
| 5 | Implementar `v_DoSVV` | Nivel 3, luego 2 |
| 6 | Barrido de SVV | Niveles 4 y 5 |
| 7 | Implementar `GetDeviatoricStrainRate` | Nivel 1 (los tres campos analíticos) |
| 8 | Implementar Smagorinsky | Nivel 3, luego 5 |
| 9 | Limitador de choque / Ducros | Nivel 6 |
| 10 | Vreman o WALE | Comparación contra Smagorinsky en 5 y 6 |

Los pasos 2 y 3 parecen triviales y son los que más gente traba. Hacelos con
los cuerpos de función **vacíos**, antes de escribir una sola línea de física.
Separar "¿está bien conectado?" de "¿está bien la fórmula?" es la diferencia
entre depurar dos cosas a la vez o de a una.
