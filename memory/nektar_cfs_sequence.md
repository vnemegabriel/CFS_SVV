# Nektar++ CompressibleFlowSolver — notes from Dev Guide ch.15, User Guide ch.9/11 (read 2026-10-04)

- Dev Guide 15.1: dU/dt = -div H, H = F(U) - G(U,gradU) (15.1); weak DG (15.2); quadrature (15.3); matrix form du/dt = M^-1[Σ B^T D_j^T Λ(wJ) H_j - B_Γ^T M_c^T Λ(w^Γ J^Γ) Ĥn] (15.4).
- RK stages (15.5-15.8); implicit: N(u)=u-s-Δt a_ii L(u)=0 (15.9), Newton (15.10), GMRES, FD Jacobian-vector (15.12), BRJ preconditioner (15.13).
- Table 15.1 functions: AdvectVolumeFlux, AdvectTraceFlux, IProductWRTDerivBase, AddTraceIntegral, MultiplyByElmtInvMass, DiffuseCalcDerivative, DiffuseVolumeFlux, DiffuseTraceFlux, AddDiffusionSymmFluxToCoeff (IP), DoOdeRhs, DoImplicitSolve, NonlinSysEvaluatorCoeff1D, MatrixMultiplyMatrixFreeCoeff, PreconCoeff.
- Fig 15.1 inheritance: EulerCFE/NavierStokesCFE -> CompressibleFlowSystem (DG calls it CompressibleFlowSolverSystem) -> AdvectionSystem -> UnsteadySystem -> EquationSystem. Members: m_advObject (+m_riemann), m_diffusion, m_artificialdiffusion, m_bndConds, m_forcing, m_intScheme.
- Table 15.2: implicit nesting Init -> time loop (UnsteadySystem) -> RK stage (TimeIntegrationScheme) -> Newton (NewtonSolver) -> residual (CFSImplicit) -> GMRES (NekLinSysIterGMRES) -> BRJ (PrecondBRJ) -> output.
- Library chapters 5-9 in Dev Guide are mostly placeholders ("xx"); Modified_A basis: vertex modes first, then bubbles (5.2.6). Basis storage: fast index = quadrature point.
- UG 9.3: Projection only DisContinuous; only WeakDG fully supported (FR only quads); DiffusionType LDGNS (penalty ∝ 1/h, LDGNSc11) or InteriorPenalty; UpwindType AUSM0-3, Average, ExactToro, HLL, HLLC, LaxFriedrichs, Roe; ShockCaptureType NonSmooth (explicit only) / Physical (NS only); ShockSensorType Modal/Dilatation; DucrosSensor; Smoothing C0.
- UG 9.4.1: ε = ε0 (h/p) λmax S (9.8); s_e = log10(<q-q̃,q-q̃>/<q,q>) (9.9); S piecewise sine (9.11), s0 = sκ - 4.25 log10 p; params Skappa, Kappa, mu0.
- UG 9.4.3: Qmin = Pexp + max(2Pexp,Pgeom)/2 + 3/2 (9.13).
- UG 11.3.5.2 (IncNS): SpectralVanishingViscosity True (exp kernel, SVVCutOffRatio, SVVDiffCoeff scaled h/p), PowerKernel, DGKernel (recommended, default coeff 1).
- Gotcha: "AdvectioType" misspelt in CylinderSubsonic_WeakDG_Implicit.xml and in UG example -> ignored.
- 1D model code/cfs1d/cfs1d.cpp: weak DG Euler, P+1 convergence; Sod P=4 Ne=40 fails without stabilisation (t=0.003), SVV c>=0.5 completes, TV ≈ 1.65 vs exact 0.875.
