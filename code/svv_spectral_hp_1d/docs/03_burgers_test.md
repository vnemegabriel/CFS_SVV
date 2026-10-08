# 3. The Burgers test of K&S Fig. 6.27

> "Solution of the inviscid Burgers equation at time $T = 0.5$ using continuous
> Galerkin (a) without, and (b) with SVV. Five equally-spaced elements spanning
> $[-1,1]$ were used, each of which contained sixteen modes. A wave cut-off of
> $M_{SVV} = 8$ and amplitude $\epsilon = 1/16$ were applied." (K&S p. 351)

Setup used here (`apps/burgers_svv.cpp`):

| item | value | note |
|---|---|---|
| domain, elements | $[-1,1]$, 5 uniform, $h = 0.4$, $J = 0.2$ | the shock sits *inside* element 2, not on an interface (the book notes an interface-aligned shock is much easier) |
| polynomial order | $P = 15$ (16 modes) | "sixteen modes" |
| quadrature | $Q = 16$ GLL (`--Q`) | classic collocation choice; over-integration available |
| initial condition | $u_0 = -\sin(\pi x)$, $u(\pm 1) = 0$ | standing shock at $x = 0$ from $t = 1/\pi$ |
| SVV | $\epsilon = 1/16$, $M_{SVV} = 8$, exponential kernel | kernel $\hat F_k = 0$ ($k \le 9$), 0.002, 0.17, 0.57, 0.85, 0.97, 1 |
| time integration | RK4, $\Delta t = 2.5 \times 10^{-3}$ (advective) | see §3.3 for the time-step dependence of the split application |

$\epsilon$ is a **physical** viscosity: the SVV term is $\epsilon\,(F * u_x)_x$ in $x$, not in $\xi$.

Diagnostics printed by the app: total variation (TV) of the sampled solution (the
exact solution has $\mathrm{TV} = 4$ for every $t$: $0 \to 1 \to -1 \to 0$; every extra unit is wiggle), $\max\lvert u\rvert$, the per-element
Legendre spectra (`*_spectrum.csv`) and, via `scripts/burgers_error.py`, the error
against the exact characteristics solution.

> **Correction (September 2026).** Earlier versions of this note were produced with a
> split application that used $c = \Delta t\,\epsilon/J$ instead of $\Delta t\,\epsilon/J^2$
> (see `docs/01` §1.4). The applied viscosity was $\epsilon J = \epsilon/5$, so the
> "Fig. 6.27(b) reproduction" of those versions was a run at $\epsilon = 1/80$, and the
> Galerkin-versus-split comparison was made at viscosities differing by a factor of 5.
> All numbers below use the corrected code. The unit test
> `split SVV step = implicit Euler of the physical element system` guards against the
> regression.

## 3.1 Verification before the shock

With the same discretisation and no SVV, error against the exact solution:

| $T$ | max error | note |
|---|---|---|
| 0.05 | $1.3 \times 10^{-10}$ | spectral accuracy |
| 0.10 | $4 \times 10^{-9}$ | |
| 0.20 | $9.2 \times 10^{-5}$ | profile steepening, shock at 0.318 |

So the spatial discretisation (basis, derivative, assembly, projection, RK4) is
correct. At $T = 0.2$ SVV with the book's parameters raises the error to
$2.9 \times 10^{-4}$ (split) and $5.2 \times 10^{-4}$ (Galerkin): the steepening profile
already has energy in modes $k > 9$, and SVV damps it.

## 3.2 Results at $T = 0.5$

`results/fig_6_27_burgers.png` shows (a) no SVV and (b) split SVV with the book's
parameters; `burgers_spectra.png` shows the element spectra;
`burgers_galerkin_vs_split.png` ($\epsilon = 1/16$) and
`burgers_galerkin_vs_split_e80.png` ($\epsilon = 1/80$) compare the two applications.

| run ($Q = 16$, conservative form) | TV | $\max\lvert u\rvert$ | comment |
|---|---|---|---|
| no SVV | 7.44 | 1.27 | wiggles around the shock, survives — Fig. 6.27(a) |
| split, $\epsilon = 1/16$ | 10.88 | 1.39 | oscillations in every element, spikes at the interfaces $x = \pm 0.2$, $\pm 0.6$ |
| Galerkin, $\epsilon = 1/16$ | 15.79 | 1.33 | same pattern, larger |
| split, $\epsilon = 1/80$ | **6.84** | 1.33 | smooth away from the shock, small oscillations at the shock: the look of Fig. 6.27(b) |
| Galerkin, $\epsilon = 1/80$ | 7.49 | 1.31 | small spikes at $x = \pm 0.2$ |
| split, $\epsilon = 1/256$ | 6.67 | 1.29 | |
| Galerkin, $\epsilon = 1/256$ | 6.55 | 1.29 | |

The sweep (`scripts/sweep_burgers.sh`) around the book's parameters:

| variant (split, $\epsilon = 1/16$ unless noted) | TV | what it teaches |
|---|---|---|
| `--Mcut 4` | 8.22 | damping more modes helps |
| `--Mcut 4 --eps 0.0125` | **6.32** | best setting found: low amplitude, wide band |
| `--eps 0.015625` ($1/64$) | 7.10 | for $\epsilon \ge 1/128$, TV grows monotonically with $\epsilon$ |
| `--form kirby` | 9.88 | eq. (6.5.17) is less dissipative, as the book says |
| `--kernel step` | 20.19 | the sharp cut-off (6.5.12) is much worse than the $C^\infty$ kernel |
| `--Mcut 12` | 11.14 | too few damped modes |
| `--eps 0.25` | 15.65 | four times the book's amplitude: worse again |
| `--Q 24 --conv nc`, no SVV → SVV | 55.5 → 13.2 | the non-conservative form is badly aliased; SVV bounds it |
| `--Q 16 --conv skew`, no SVV → SVV | 31.1 → 13.1 | same story |
| `--Q 24`, no SVV → SVV | 22.8 → 10.1 | over-integrating the conservative form is *worse* than $Q = 16$ here |
| `--Q 17`, no SVV → SVV | blow-up → 17.0 | an odd quadrature count is the least stable option; SVV keeps it bounded |
| Galerkin, `--Q 24` | 13.2 | |
| Galerkin, `--eps 0.25` | 24.5 | the Galerkin artefact grows with $\epsilon$ faster than the split one |

## 3.3 Time-step dependence of the split application

The split application is a first-order (Lie) splitting: RK4 for convection, one
implicit-Euler SVV step per time step. Its result depends on $\Delta t$; the Galerkin
application does not (once $\Delta t$ is below its stability limit).

| $\Delta t$ | split, $\epsilon = 1/16$ | Galerkin, $\epsilon = 1/16$ | split, $\epsilon = 1/80$ | Galerkin, $\epsilon = 1/80$ |
|---|---|---|---|---|
| $2.5 \times 10^{-3}$ | 10.91 | unstable | 6.84 | 7.49 |
| $1.25 \times 10^{-3}$ | 12.38 | 15.60 | | |
| $6.25 \times 10^{-4}$ | 13.28 | 15.79 | 7.02 | 7.49 |
| $3.125 \times 10^{-4}$ | 13.77 | 15.79 | | |
| $1.5625 \times 10^{-4}$ | 14.03 | 15.79 | 7.07 | 7.49 |

The split TV converges at first order to $\approx 14.3$ ($\epsilon = 1/16$) and
$\approx 7.1$ ($\epsilon = 1/80$). At the advective time step, part of the split's
advantage is splitting error. The converged split and Galerkin results still differ,
because the two applications are different semi-discretisations: the split one
freezes the vertices during the SVV step, the Galerkin one lets the consistent mass
matrix move them.

## 3.4 What the test shows

1. **The operator is correct; the book's amplitude does not reproduce Fig. 6.27(b)
   in this code.** With $\epsilon = 1/16$ as a physical viscosity, both applications
   raise TV above the no-SVV value, and the spectra plot shows modes 11–13 *growing*
   by about two orders of magnitude in elements 1 and 3 (the neighbours of the shock
   element). The figure's look is obtained at $\epsilon \approx 1/80$. Either the
   book's $\epsilon$ is defined on the reference element (then $\epsilon_{phys} =
   \epsilon J$ or $\epsilon J^2$, depending on the convention), or Kirby's code differs
   in another respect. The text does not settle it.
2. **The interface artefact is not specific to the Galerkin application.** At
   $\epsilon = 1/16$ the split application also produces spikes at the element
   interfaces and high-mode growth in the neighbouring elements. The mass-matrix
   mechanism of `docs/01` §1.4 explains why Galerkin is worse at equal $\epsilon$ and
   equal $\Delta t$ convergence (15.8 against 14.3), not the artefact itself. Open
   question: what injects high-degree content in elements that SVV does not couple?
3. **SVV bounds aliasing-unstable discretisations** (nc, skew, $Q = 17$, $Q = 24$),
   but at the book's amplitude it leaves them worse than the stable $Q = 16$
   conservative form without SVV.
4. **Low amplitude and wide band work best here**: TV grows monotonically with $\epsilon$
   from $1/128$ to $1/4$ (split: 6.63 at $1/128$, 6.67 at $1/256$), and `--Mcut 4 --eps 0.0125` gives the lowest TV of the sweep (6.32).

## 3.5 Exercises

* Change `--Mcut` from 2 to 14 at `--eps 0.0125` and plot TV; find the value that
  minimises it. Compare with the $5\sqrt{P} \approx 19$ rule of Maday et al.
  (impossible here: $P = 15$).
* Run `--nel 10 --P 7` (same number of dofs). Without SVV this blows up at
  $t \approx 0.41$; with `--Mcut 4` it survives (TV 7.14). Why is the
  low-order/many-element version less robust without SVV?
* Run a single-domain Legendre case, the setting of the SVV theory:
  `--nel 1 --P 79 --Q 81 --Mcut 45 --eps 0.0125`. Compare with `--eps 0.05`.
* Investigate open question 2: run split with `--eps 0.0625 --dt 1.5625e-4` and plot
  the spectra of elements 1 and 3 over time. Does the growth start at the shock
  formation time $t = 1/\pi$?
* Replace implicit Euler in `SVVBubbleFilter::apply` by Crank–Nicolson and repeat the
  table of §3.3. Does the split result still depend on $\Delta t$ at first order?
