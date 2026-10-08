> **SIEMPRE PENSAR: que es importante y clave de esto que leo?**


## 1/8 
Trabajando en el workflow para poder correr casos de forma más sencilla

## 31/8
- CASOS
	- ver como plottear conjunto de checkpoints de forma fluida para ver la simulación completa
	- seteo de configuración supersónica
		- como simplifico la configuración??
	- 
- 

## 25/8

Pasando en limpio **notas escritas**
- ==ESTRUCTURAR== y ==SIMPLIFICAR== casos como corresponde 
	- 3 casos de paper de Nektar `(code/nektar_cases/nektar_paper)`
		- cylinder
		- naca0012
		- t106c (turbina)
	- 3 casos de paper de SVV de Mateo Gabin `(code/nektar_cases/svv_paper)`
		- shu-osher
		- tgv
		- escalon
	- Aconcagua 
		- Malla en `meshingMaster/aconcagua`
			- Plantear convergencia de la malla (no irse de mambo)
--- 
## 14/7

Reunion con Pablo
- Hacer clases de 20 min para los martes, resumiendo y explicando lo que se hizo en la semana
- Revisar plan de tesis, que voy a modificar?
- ordenar bibliografia y resultados que buscamos, no estar tanteando
- arrancar por lo basico, navier stokes. cual es el problema en compresible?
	- como se resuelve el problema de la densidad?
	- que modelos se usan?
	- por que tu alternativa?
	- que compone tu alternativa?

Leyendo Nodal DG methods (libro)

cap 2.

seccion 2.1 

space filling triangulation -> Operacion matematica de creacion geometrica de elementos no superpuestos que sirve como base para definir una malla

# 20/7

cap 3 (voy mechando cn cap 2)

## 21/7 

Explicación rápida de: Por qué usamos una matriz Vandermonde?
[[5b3c34a6e02d7eab925129a30e7ed4e0_MD5.jpg|Open: Pasted image 20260721193724.png]]
![[5b3c34a6e02d7eab925129a30e7ed4e0_MD5.jpg]]
Mostrar visualizacion de ejemplo 3.1
## 22/7

Arrancar desde el principio.

- implementacion 1d de ns con dgsem
- modificacion plan de tesis 

## 25/7

Lectura y avance de NDG hasta generacion de malla y de mapeo

## 1/8

Revisar : [High Order Positivity-Preserving Entropy Stable Discontinuous Galerkin Discretizations](https://www.youtube.com/watch?v=x33D171lXjc)
![[b8905395c56d59711f7bf31f0ef61334_MD5.jpg]]

Deje el video en minuto 20. Terminar para entender bien la teoría
#### Proximos pasos

Reescribir plan de tesis a **mano** para tener claro en la cabeza:
- Que son los metodos espectrales
- Por que los uso
- Que aportan en este caso
- Que es SVV
- Que es Nektar++
- Que voy a hacer de aca a diciembre
- Que casos voy a probar
- Que comparaciones practicas contra FVM voy a hacer (OpenFOAM)

Armar to-do list en partecitas logrables a diario

## 10/8

### Que voy a hacer?

Agregar el filtro SVV a la biblioteca de funciones de CompressibleFlowSolver. 

### Que busco probar? 

Que para un alto rango de


