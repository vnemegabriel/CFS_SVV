# Cadena: formación del choque → espacio espectral → discretización → SVV → calibración (2026-10-04)

Marcas: [V-C] texto leído por Claude; [V-bib] registro bibliográfico confirmado (Crossref/arXiv), texto no leído; [F] no verificado; [X] error detectado.
Pendiente: K–S cap. 10 completo (el PDF de 34 MB no pasa por el puente; sec. 10.1 ya está en karniadakis_10_1.md).

## 1. Formación analítica del choque en aire viscoso

- Empinamiento inviscido (Burgers): $u_t+uu_x=0$, a lo largo de características $u_x=\dfrac{u_0'(\xi)}{1+u_0'(\xi)\,t}$, explota en $t^*=-1/\min_\xi u_0'(\xi)$. Sistemas 2×2: Lax (1964); n×n genuinamente no lineales: John (1974); 3D: Sideris (1985), Christodoulou (2007), Buckmaster–Shkoller–Vicol (2023).
- Onda acústica débil en gas ideal: velocidad característica $c_0+\beta u'$, $\beta=(\gamma+1)/2=1.2$ en aire.
- Con viscosidad y conducción (reducción de Lighthill 1956 [F]): $u_t+\beta u u_x=\tfrac{\delta}{2}u_{xx}$, $\delta=\nu\left(\tfrac43+\tfrac{\mu_B}{\mu}+\tfrac{\gamma-1}{\Pr}\right)$ [F, forma estándar, verificar contra Lighthill].
- Perfil de Taylor (1910): $u=\bar u-\tfrac{\Delta u}{2}\tanh\!\left(\tfrac{\beta\Delta u\,x}{2\delta}\right)$, espesor de pendiente máxima $L=4\delta/(\beta\Delta u)$ [derivación propia desde Burgers].
- Choque fuerte: perfil de Becker (1922) exacto con $\Pr=3/4$ (entalpía total constante a través del choque [V-C, Persson–Peraire sec. III.B.2]); existencia general: Gilbarg (1951) [V-bib, alcance de hipótesis F].
- Escala de espesor (Persson–Peraire ec. 19, citando Liepmann–Roshko) [V-C]: $\dfrac{\delta_s\,\Delta u\,\rho^*}{\mu^*}\approx 1$, $\Delta u=\dfrac{2u_\infty(M_\infty^2-1)}{(\gamma+1)M_\infty^2}$.
- Saltos de Rankine–Hugoniot independientes de $\mu,\kappa$ (integrar la forma conservativa entre estados uniformes). Límite $\mu\to0$: solución de entropía (Hoff–Liu 1989 isentrópico; Bianchini–Bressan 2005 viscosidad artificial uniforme).
- Espesor físico del orden de algunos recorridos libres medios [V-C, Pirozzoli 2011 p.165]; Navier–Stokes ya difiere del experimento a $M=1.55$ [V-bib, resumen de Alsmeyer 1976].
- Consecuencia: el espesor numérico del choque lo fija la estabilización, nunca $\mu$ físico.

## 2. Traducción al espacio espectral

- Galerkin–Fourier de Burgers [V-C, Tadmor ec. 5.1]: $\dfrac{d\hat u_k}{dt}+\dfrac{ik}{2}\sum_{p+q=k}\hat u_p\hat u_q=0$. Empinamiento = transferencia por tríadas hacia $k$ alto.
- Sin disipación se conserva $\int u_N^2$ (Tadmor ec. 1.9); la solución de entropía disipa $\int u^2$ en el choque ⇒ el método espectral puro no puede converger a ella (Tadmor sec. 1; Theorem A1: tampoco con post-proceso por convolución) [V-C]. Legendre: mismo fracaso (Maday et al. Fig. 6.1) [V-C].
- Decaimiento de coeficientes (base ortonormal) [derivación; P–P sec. II para el caso continuo]: salto $|\hat u_k|\sim k^{-1}$; continuo con quiebre $\sim k^{-2}$; analítico, exponencial. Sensor $S_e=\|u-\hat u\|^2/\|u\|^2$: salto $\sim p^{-2}$, quiebre $\sim p^{-4}$ (P–P) [V-C].
- Perfil $\tanh(x/\ell)$: transformada $\propto \ell\,\mathrm{csch}(\pi k\ell/2)$, decae como $e^{-\pi k\ell/2}$; el perfil está resuelto si $\ell\gtrsim h/p$ (resolución $\delta\sim h/p$, P–P sec. I [V-C]).

## 3. Problemas al discretizar

- Gibbs [V-C, Gottlieb–Shu sec. 1]: error $O(1/N)$ lejos del salto; sobrepico que no decae con $N$.
- Monotonicidad: Godunov (1959) [F]; esquemas monótonos a lo sumo de primer orden (Harten et al. 1986, vía Pirozzoli p.170 [V-C]). Alto orden sin oscilación exige disipación no lineal (dependiente de la solución) ⇒ sensor.
- Aliasing de términos no lineales (Kirby–Sherwin sec. 1) [V-C]; positividad: la estabilidad entrópica requiere $\rho>0$ (Mateo-Gabín sec. 1) [V-C].

## 4. SVV

- Tadmor (1989) [V-C]: $\partial_t u_N+\partial_x P_N(u_N^2/2)=\varepsilon\,\partial_x\big((I-P_m)\partial_x u_N\big)$; $\varepsilon\sim N^{-2\beta}$, $m\sim N^{\beta}$, $0<\beta<1/4$ ⇒ convergencia a la solución de entropía (Thm 4.1), condicionada a cota $L^\infty$. Sistemas: viscosidad aplicada a variables de entropía $v=\partial U/\partial u$ (sec. 5). Experimento: $\varepsilon m=0.25$ (lectura OCR).
- Maday–Ould Kaber–Tadmor (1993) [V-C]: Legendre, Thm 5.2 (exponentes con OCR ambiguo); práctica $\varepsilon_N=N^{-1}$, $m_N=5\sqrt N$, $\hat Q_k=\exp\!\big(-(k-N)^2/(k-m_N)^2\big)$. La solución SV converge en orden bajo; la precisión espectral se recupera sólo con post-proceso.
- Kirby–Sherwin (2006) [V-C]: base modal ortogonal necesaria para operador simétrico semidefinido positivo; filtro ec. 15 por orden total; Burgers $P=15$, $P_{cut}=7$, $\varepsilon=1/16$. [X] Atribuyen a Tadmor "$P_{cut}\approx0.25$"; en Tadmor el 0.25 es el producto $\varepsilon m$ (Mateo-Gabín ec. 66 lo escribe $\varepsilon M=0.25$).
- Mateo-Gabín et al. (2022) [V-C]: núcleo de potencia $\hat F_i=(i/N)^{P_{SVV}}$ (Moura); "SVV filtering is not enough to control the oscillations near the strong shock" → sensor de gradiente de densidad; FFS Mach 3, $N=7$: $P_{SVV}=4$ fuera del choque, $0$ en el choque, $\mu=5\times10^{-4}$. App. A: $P_{SVV}=O(1)$ supersónico, $O(10^{-2})$ turbulento; sin estrategia de ajuste, prueba y error. Archivo local mal nombrado (Gabin2008 → 2022).
- Nektar++ (UG 9.4.1): $\varepsilon=\varepsilon_0(h/p)\lambda_{max}S$, $s_0=s_\kappa-4.25\log_{10}p$.

## 5. Calibración

- Salto y posición: forma conservativa + convergencia (Lax–Wendroff) + selección entrópica. La prueba de Tadmor para sistemas usa variables de entropía; con SVV sobre variables conservativas esa prueba no aplica [decidir].
- Espesor objetivo $\ell=c\,h/p$; amplitud por ec. 19 de P–P: $\mu^*\approx\rho^*\Delta u\,c\,h/p$.
- Comparación con la escala de Nektar++ ($\lambda_{max}=u_\infty(1+1/M)$): $\Delta u/\lambda_{max}=\dfrac{2(M-1)}{(\gamma+1)M}$ → 0.14 ($M=1.2$), 0.42 ($M=2$), 0.56 ($M=3$). Con $\varepsilon_0$ fijo, el choque débil recibe ~4× más viscosidad relativa que el de $M=3$; coherente con el ajuste por Mach que reporta P–P [derivación propia, verificar numéricamente].
- Las leyes asintóticas no aplican por elemento a $p=3$–$8$: $m_N=5\sqrt p\ge p$ para $p\le25$. Régimen pre-asintótico ⇒ calibración empírica de núcleo, $P_{cut}$ y amplitud.
- Umbral del sensor: separar salto ($\log_{10}S\approx-2\log_{10}p$) de suave ($\le-4\log_{10}p$).
- Verificación: perfil viscoso de Becker ($\Pr=3/4$, $\mu$ suficiente para resolver $\delta_s$) contra solución exacta; espesor numérico vs $h/p$ y $M$; Sod, Shu–Osher, FFS $M=3$.

## 6. Herramientas faltantes

- Estabilidad entrópica sin SBP/split-form (Nektar++ no lo tiene): proyección a variables de entropía [F].
- Limitador de positividad: Zhang–Shu, JCP 2010, doi 10.1016/j.jcp.2010.08.016 [V-bib].
- Subcelda FV (Hennemann et al. 2021, JCP 426:109935) [V-bib]: requiere DGSEM con SBP.
- Discriminar choque/contacto: Ducros/dilatación (existen en Nektar++, sólo con ShockCaptureType Physical).
- Recuperar precisión espectral: Gegenbauer / filtros unilaterales (Gottlieb–Shu sec. 4; Maday sec. 7) [V-C].
- Elección de núcleo por análisis de autosoluciones: Moura et al. 2016 [V-bib, no está en raw/].
