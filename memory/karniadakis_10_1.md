# Karniadakis & Sherwin, 2nd ed., §10.1 Conservative formulation (pp. 515–524, PDF pp. 534–543)

- (10.1.1) u_t + f(u)_x = 0.
- (10.1.2) Cell average: ū_i = 1/(x_{i+}-x_{i-}) ∫ u dx.
- (10.1.3a) dū_i/dt + [f(u_{i+}) - f(u_{i-})]/Δx_i = 0; fluxes evaluated with de-averaged (reconstructed) point values at cell ends. Staggered grid.
- (10.1.4) Fourier: ā_k = σ_k a_k, σ_k = sin(kΔx/2)/(kΔx/2), Δx = π/P.
- (10.1.5) Chebyshev: cell centres = Gauss points x_i = cos((i-1/2)Δθ), 1≤i≤P; cell ends = Gauss–Lobatto x_{i±} = cos(i±Δθ), 0≤i±≤P; Δθ = π/P.
- (10.1.6) T̄_0 = 1, T̄_1 = σ_1 U_1/2, T̄_k = [σ_k U_k - σ_{k-2} U_{k-2}]/2; σ_k = sin((k+1)Δθ/2)/((k+1) sin(Δθ/2)); U_k = T'_{k+1}/(k+1).
- (10.1.7) Spectral element: u^e = Σ u_n h_n, h_n = (2/P) Σ_p T_p(x_n)T_p(x)/(c̄_n c̄_p), c̄ = 2 at n=0,P.
- (10.1.8a) ū = A u, A_in = h̄_n(x_i). (10.1.8b-c) Gauss–Chebyshev Lagrange g_j = T_P/(T'_P(x_j)(x-x_j)).
- (10.1.10–13) Reconstruction: degree P-1 polynomial from P averages via de-averaged interpolants G_j (λ_p^j coefficients in U_p basis).
- (10.1.14) Interfacial constraint: one extra value per element fixed by upwinding, u_γ from upwind element; Burgers: V = (u^e_P + u^{e+1}_0)/2.
- §10.1.4 Non-oscillatory: Algorithms R (reconstruction) and A (averaging) split a single jump into step + smooth part (refs [83], [438]); test fn (10.1.15), jump at x_s = 5, K = 5, P = 20/40/80, Vandeven filter.
- §10.2 Monotonicity: FCT (10.2.1), local projection limiting (10.2.2).

Implementation: code/karniadakis/hyperbolic_10_1.cpp — single Chebyshev domain, exact inverse of [inflow row; A] for reconstruction, RK4. Spatial error spectral (≈1e-11 at P=24, linear advection).

## Prerequisites for §10.1 (read 2026-10-04)
- §2.3.1 (2.3.1): affine map x = (1-ξ)/2 x_{e-1} + (1+ξ)/2 x_e, J = h_e/2; averages invariant.
- §2.3.4.1 (2.3.9–10): h_p(x_q)=δ_pq, h_p = g/(g'(x_p)(x-x_p)); GLC g=(1-ξ²)T_P', Gauss g=T_P. App A (A.1.11): Chebyshev sum form.
- App A.1: T_n = 2^{2n}(n!)²/(2n)! P_n^{-1/2,-1/2}; T_n'=nU_{n-1}; T_n=(U_n-U_{n-2})/2; ∫T_n = ½[T_{n+1}/(n+1) - T_{n-1}/(n-1)].
- §2.4.1, App B: Gauss exact to 2Q-1, Gauss–Lobatto 2Q-3; Chebyshev points explicit (B.2). §2.4.1.2 Table 2.2: Q_min = P+2 (linear), 3P/2+2 (quadratic), 2P+2 (cubic).
- §4.1.6 (4.1.38): sum factorisation O(P^4)→O(P^3) in 2D.
- §6.3.1 (6.3.2): Δt ≤ α_im/(C P²), α_im=0.723 AB3; c_λ≈0.2; §6.3.2.1 Chebyshev Dirichlet max|λ|≈0.089P².
- §6.5: Godunov; §6.5.1 Vandeven filter (6.5.4); §6.5.2 SVV (Tadmor, Burgers 6.5.8); App E: A = R D L characteristic upwinding.
- Measured: staggered operator (linear advection, inflow left) max|λ| ≈ 0.22 P², all Re λ < 0; RK4 limit C_max ≈ 2.6 w.r.t. Δξ_min.
- Worked P=2: A = [5 8 -1; -1 8 5]/12, R = [1 0 0; -1/2 5/4 1/4; 1 -2 2].
