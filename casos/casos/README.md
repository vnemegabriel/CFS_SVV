# Casos

Corre matrices de simulaciones, vigila cada corrida, clasifica sus fallas sin cortar la caso
y rehace los resultados normalizados. Python 3 con la biblioteca estándar; matplotlib es opcional
y sólo agrega figuras.

Casos, en el orden de niveles de `TESIS/CASOS.md`:

| Casos               | Nivel              | Resultado                                                |
| ------------------- | ------------------ | -------------------------------------------------------- |
| `verificacion/`     | V2–V4              | Normas de error y criterio binario por caso              |
| `calibracion_tubo/` | C2                 | Normas de error de las sesiones del tubo de choque       |
| `naca0012/`         | Validación y costo | Curva costo–precisión (`TESIS/CURVA_COSTO_PRECISION.md`) |
| `aconcagua/`        | Aplicación         | Cd contra Mach y contra telemetría                       |

Las plantillas de `verificacion/` y `calibracion_tubo/` salen del test suite de Nektar++ con su
`generar.sh`, que deja a la vista qué se cambió de cada sesión original.

## Archivos

| Archivo | Hace una cosa |
|---|---|
| `casos.py` | Línea de comandos |
| `cfg.py` | Lee `maquinas.ini`, `case.ini`, `matrix.txt`; resuelve los valores de una corrida |
| `runner.py` | Una corrida: tomar, preparar, lanzar, vigilar, registrar |
| `parsers.py` | Salidas de cada solver: fuerzas, tiempos, celdas, firmas de falla |
| `post.py` | Historial → punto de parada → `result.ini` normalizado |
| `report.py` | `resultados/` de la caso, rehechos desde cero |
| `num.py` | Interpolación, criterio de parada, PCHIP, Richardson |
| `telemetry.py` | Cd de vuelo por desaceleración |
| `render.py` | Plantillas con `@nombre@` |
| `maquinas.ini` | Una sección por máquina |
| `<caso>/case.ini` | Qué correr y cómo post-procesarlo |
| `<caso>/matrix.txt` | Una corrida por línea |
| `<caso>/resultados/` | `results.csv`, `resumen.txt` y las tablas de la caso |

## Uso

```bash
python3 casos.py plan   naca0012        # matrix.txt desde [sweep]; revisarla a mano
python3 casos.py check  naca0012        # verifica todo sin correr nada
python3 casos.py run    naca0012        # en serie, en esta máquina
python3 casos.py status naca0012
python3 casos.py post   naca0012        # rehace resultados/
```

`run --solo 'nek_dg_M0.8_*'` restringe por patrón. `run --retry` repite las que no terminaron `ok`.
Cortar con Ctrl+C marca la corrida en curso como `interrumpida` y mata el solver.

En un cluster:

```bash
python3 casos.py submit naca0012        # un trabajo por corrida
python3 casos.py post   naca0012 --cada 600
```

## Qué deja cada corrida

`<runs_root>/<caso>/<id>/`

| Archivo | Contenido |
|---|---|
| `status` | Una línea: `estado \| fecha \| host \| motivo` |
| `run.ini` | Valores y comandos exactos con que se lanzó |
| `log` | Salida del solver |
| `prepare.log` | Salida de la preparación |
| `case/` | Directorio de trabajo del solver |
| `case.anterior/` | El `case/` de la corrida previa con el mismo id. Se guarda uno solo |
| `result.ini` | Resultado normalizado. Existe también en corridas fallidas, con lo que haya |
| `runner.err` | Traza, sólo si falló el propio programa |
| `lock` | Mientras corre. Su fecha de modificación es el latido |

## Estados

| Estado | Significa |
|---|---|
| `pendiente` | Nunca se lanzó |
| `encolada` | Enviada a la cola |
| `corriendo` | Con latido reciente |
| `corriendo?` | Sin latido: el runner murió. `casos.py unlock` la libera |
| `ok` | Terminó y cumplió el criterio de parada |
| `sin_converger` | Terminó sin cumplirlo |
| `divergio` | NaN en el log o en las fuerzas |
| `fallo` | Error de configuración, de preparación, del solver o código de salida distinto de cero |
| `vencida` | Superó `timeout_s` |
| `estancada` | Sin salida nueva durante `stall_s` |
| `interrumpida` | Ctrl+C, señal de la cola, o `unlock` |

Un estado terminal no se repite solo: hace falta `--retry`. Cada repetición arranca de cero.

## Normalización

Todo resultado sale en las mismas unidades, sin importar solver ni máquina.

| Magnitud | Definición |
|---|---|
| Cd, Cl | Fuerza total del solver proyectada con `alpha` y `beta`, dividida por `0.5·rho·U²·aref` de `[values]`. No se usa la normalización propia del solver |
| Tiempo | Unidades convectivas `t·U/lref`, o iteraciones si `unit = iter` |
| Parada | Primer instante en que Cd y Cl varían menos que `tol_cd` y `tol_cl` durante la última `window` |
| Resultado | Promedio de Cd y Cl en esa ventana |
| Segundos de solver | Hasta la parada, sin escritura a disco. Nektar++: bloques `Steps` menos `Writing`. OpenFOAM: `ClockTime` desde el inicio del lazo, resolución de 1 s |
| Horas-núcleo | procesos × segundos de solver / 3600 |
| Costo normalizado | procesos × segundos de solver / `taubench_s` de la máquina |
| Grados de libertad | Expresión `dof` de `[post]` con `celdas` contadas de la malla |
| Tamaño efectivo | grados de libertad^(−1/`dim`) |

`taubench_s` es el tiempo de TauBench en esa máquina. Mientras falte, el costo normalizado queda
vacío y `resumen.txt` lo avisa; al completarlo, `casos.py post --reprocesar` lo llena.

## `case.ini`

| Sección | Claves |
|---|---|
| `[case]` | `solver` (`nektar` u `openfoam`), `id` con `{columnas}` |
| `[sweep]` | Una clave por columna; `plan` arma todas las combinaciones. `mach,alt = 0.8,3000 1.8,1200` agrupa columnas que van juntas. Un valor `PENDIENTE` frena `check` |
| `[values]` | Expresiones Python en orden; ven las columnas, las anteriores y `var_*` de la máquina. Obligatorias: `rho`, `U`, `aref`, `lref` |
| `[run]` | `prepare`, `command`, `watch`, `probe`, `timeout_s`, `stall_s`, `poll_s` |
| `[post]` | `kind` (`forces` o `errores`), `forces`, `log`, `mesh`, `dof`, `dim`, `unit`, `window`, `tol_cd`, `tol_cl` |
| `[watch]` | `divergio`, `fallo`: expresiones regulares extra |
| `[report]` | `kind` (`cost_curve`, `cd_mach` o `errores`) y sus claves |
| `[criteria]` | Sólo con `errores`, una línea por caso: `linf <= X`, `igual stab=dg tol=X`, `orden h variable=V margen=X` |
| `[telemetry]` | Mapa de columnas del archivo de vuelo, masa, área. Con `col_accel` el Cd sale del acelerómetro; sin él, de derivar la velocidad |

En `[run]` y `[post]`, `{nombre}` se reemplaza por el valor de la corrida, más `{run}`, `{case}`,
`{case}` y `{tools}`. Un valor que no se puede calcular hace fallar `check`, no la corrida.

## Varias máquinas

Cada máquina tiene su sección en `maquinas.ini`. Para juntar resultados alcanza con copiar tres
archivos por corrida:

```bash
rsync -a --include='*/' --include='status' --include='run.ini' --include='result.ini' \
      --exclude='*' cluster:/scratch/USUARIO/runs/ ~/runs_cluster/
python3 casos.py post naca0012 --desde ~/runs ~/runs_cluster
```

Sin archivos de fuerzas, `post` usa el `result.ini` copiado tal cual.

## Pruebas

```bash
python3 -m unittest discover -s tests -v
```

Parsers contra salidas reales de Nektar++ y OpenFOAM (`tests/fixtures`, incluidas una divergencia
de cada uno), un solver falso que falla de cada forma posible, interrupción con Ctrl+C,
Richardson, PCHIP, y Cd recuperado de un vuelo sintético.
