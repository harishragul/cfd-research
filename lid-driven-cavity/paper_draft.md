# Lid-Driven Cavity Flow: A Comparative Study of Projection-Method Discretisations and a Physics-Informed Neural Network, Benchmarked Against Ghia et al. (1982)

## Introduction

In my earlier note (Rajaramaduraikarthik, 2026) I built a custom Python solver for the 2D incompressible Navier-Stokes equations using the projection method, and validated it against the classic lid-driven cavity benchmark of Ghia et al. (1982). That solver used a collocated grid with first-order upwind convection, and it agreed well with the benchmark except for a small overshoot in the v-velocity trough near x = 0.8, which I attributed to numerical diffusion from the upwind scheme. In the Discussion I noted that replacing upwind with a higher-order TVD scheme should reduce this diffusion, but that doing so on the same collocated grid caused the pressure Poisson solve to misbehave — likely, I thought, an interaction between the flux limiter and the checkerboard pressure modes that collocated grids are known to suffer from (Rhie & Chow, 1983). I said a staggered MAC grid (Harlow & Welch, 1965) with a direct sparse solver would be the more robust way forward, and that this investigation was ongoing.

This note is that investigation, finished. Before touching the 2D solver again, I went back to the 1D advection groundwork mentioned in the original note and ran a proper screening study: five convection schemes (upwind and four TVD flux limiters) on a simple 1D advection test, so I could pick a limiter for the 2D work based on measured behaviour rather than just picking Van Leer again out of habit. I then built three more 2D solver variants on top of the original: the TVD scheme on the original collocated grid, the original upwind scheme on a staggered MAC grid, and TVD convection on a staggered MAC grid — four solvers in total, spanning every combination of {grid arrangement} × {convection scheme}. I also wanted to see how a completely different approach compares on the same problem, so I built a fifth solver: a Physics-Informed Neural Network (PINN, Raissi et al., 2019) that solves the same steady Navier-Stokes system by minimising a physics-residual loss instead of marching a grid forward in time. All five 2D solvers are validated against the same Ghia et al. (1982) tabulated data at Re = 100, using the same quantitative error metric, so the results are directly comparable.

## Problem Description

The lid-driven cavity problem is a 2D incompressible viscous flow problem where a unit square cavity contains fluid inside with three fixed walls and one moving wall on top (like the top wall moving like a treadmill conveyor belt). The moving top wall drags the fluid near it, and because the flow has nowhere to go, it recirculates and forms a primary vortex. It remains the standard test case for incompressible viscous flow solvers for the same reason I picked it the first time: the geometry is simple, but the flow contains recirculation, pressure gradients, and boundary layer effects that genuinely exercise a solver.

## Numerical Methods

All four classical solvers use the same projection (fractional-step) method described in the original note: a momentum predictor that ignores pressure and produces an intermediate velocity **u**\*, a Poisson solve for the pressure correction that makes the velocity field divergence-free, and a corrector step. What changes between the four variants is the grid arrangement and the convection discretisation.

**Grid arrangement.** The original solver stores u, v, and p all at the same 129×129 cell-centre locations (a *collocated* grid). Collocated grids are simple to implement but are known to admit spurious checkerboard pressure oscillations, because the discrete pressure gradient at a node doesn't "see" the pressure at that same node — a decoupling first analysed by Rhie & Chow (1983). The staggered variants use a Marker-and-Cell (MAC) grid (Harlow & Welch, 1965) instead: pressure lives at 129×129 cell centres, u lives on the vertical cell faces, and v lives on the horizontal cell faces. This staggering removes the checkerboard mode by construction, at the cost of needing interpolation whenever a variable is needed at a location it isn't natively stored.

**Convection scheme.** The original solver discretises the nonlinear convection term with first-order upwind differencing, which is unconditionally stable but numerically diffusive — it smears sharp gradients. The TVD variants replace this with the Van Leer flux limiter (van Leer, 1974), which blends toward a second-order central scheme in smooth regions of the flow while falling back to first-order upwind near steep gradients, avoiding the oscillations a naive second-order scheme would introduce.

**Pressure solve.** The two collocated solvers use Jacobi iteration on the pressure Poisson equation (max 20,000 iterations per time step, tolerance 1e-6, Neumann boundary conditions, pressure pinned to zero at one interior point). The two staggered solvers instead assemble the discrete Poisson operator as a sparse matrix once and LU-factorise it up front (`scipy.sparse.linalg.factorized`), so each time step's pressure solve is a single fast back-substitution rather than an iterative loop — this is the "direct sparse pressure solver" I said the staggered approach would let me use.

**PINN.** The fifth solver replaces the grid entirely. A fully-connected neural network (6 hidden layers, originally 64 units wide, widened to 100 in this note's final round — see Results) maps a point (x, y) directly to (u, v, p). Instead of marching a grid forward in time, the network is trained by gradient descent (Adam, then L-BFGS) to minimise a composite loss: the residual of the steady Navier-Stokes equations at ~20,000 randomly-sampled interior collocation points, plus a penalty for violating the no-slip/moving-lid boundary conditions, plus a small term pinning the pressure gauge. All spatial derivatives needed for the physics residual are obtained by automatic differentiation through the network (`torch.autograd.grad` with `create_graph=True`), not by finite differences.

## Choice of Flux Limiter: A 1D Screening Study

Van Leer wasn't the only limiter I could have picked, and I didn't want to reuse it in the 2D solver just out of habit. Before touching the 2D code, I went back to the 1D linear advection groundwork mentioned in the Introduction and ran a proper comparison: five convection schemes — first-order upwind, and four TVD flux limiters (minmod, Van Leer, superbee, and MC) in the standard Sweby (1984) flux-limiter framework — advecting a square pulse at constant velocity c = 1 on a periodic unit-length domain (N = 100, CFL = 0.5) for exactly one full domain transit. Because the domain is periodic and the pulse travels exactly one domain length in that time, the exact solution at the final time is just the initial condition again, which makes the L2 error trivial to compute without needing an analytical solution to the advected profile itself.

**Table 1. 1D advection: L2 error and monotonicity, one full periodic transit**

| Scheme | L2 error | Overshoot (max − 1) | Undershoot (min) |
|---|---|---|---|
| Upwind | 0.1817 | −0.0341 (none) | 0.0000 (none) |
| Minmod | 0.1365 | 0.0574 | −0.0556 |
| Van Leer | 0.1317 | 0.0863 | −0.0817 |
| MC | 0.1303 | 0.1073 | −0.1017 |
| Superbee | 0.1251 | 0.1507 | −0.1466 |

The results show the classic accuracy/monotonicity trade-off these limiters are known for. Upwind has by far the worst L2 error — about 38% higher than any of the TVD schemes — because it smears the pulse's sharp edges, but it is perfectly monotonic (no overshoot or undershoot at all) by construction. Superbee sits at the opposite extreme: it has the lowest L2 error of the five (sharpest resolution of the pulse), but also by far the worst overshoot and undershoot — over twice Van Leer's. This is a well-documented property of superbee: it lies at the most "compressive" edge of the Sweby TVD region and can over-sharpen gradients into new local extrema even while technically remaining TVD in the strict sense.

Van Leer is not the best on either single axis — MC and superbee both beat it on L2 error, and minmod beats it on monotonicity — but it gives up very little on either: about 28% lower L2 error than upwind, while keeping overshoot/undershoot noticeably smaller than MC or superbee. Given that the 2D lid-driven cavity solver already has a known sensitivity around the collocated-grid pressure solve (see Discussion), I judged the extra sharpness superbee or MC would buy in 1D wasn't worth the added oscillation risk in 2D, and kept Van Leer for all four TVD/staggered variants below.

## Grid Setup and Discretisation

All five solvers use a unit square domain (L = H = 1) at Re = 100. The four classical solvers use a 129×129 grid (either 129×129 collocated nodes, or 129×129 pressure cells with 129×130 and 130×129 velocity faces on the staggered grid) and explicit forward-Euler time integration, with the time step chosen from the combined CFL/diffusion stability limit — the same Δt = 0.8/(U/Δx + U/Δy + 4ν/Δx²) form as the original note, with a tighter safety factor (0.4 instead of 0.8) for the TVD variants, which need CFL ≤ 0.5 per direction for the limiter to remain stable. The PINN uses no grid at all — collocation points are drawn from a 2D Halton low-discrepancy sequence for good space-filling, with roughly a fifth of the points deliberately concentrated near the moving lid (y → 1), where the velocity boundary layer is sharpest.

## Results

Table 2 reports the maximum absolute error between each solver's centreline velocity profiles and the Ghia et al. (1982) tabulated values at Re = 100, using the same methodology for every solver: the simulated u(x=0.5, y) and v(x, y=0.5) profiles are cubic-interpolated onto the exact Ghia sample points, and the largest absolute difference is taken. This is a stricter and more precise check than the purely visual comparison in the original note.

**Table 2. Maximum error against Ghia et al. (1982), Re = 100**

| Solver | Grid | Convection | max\|u − u_Ghia\| | max\|v − v_Ghia\| |
|---|---|---|---|---|
| Projection method (original note) | Collocated | Upwind | 0.0103 | 0.0091 |
| Projection method | Staggered MAC | Upwind | 0.0107 | 0.0046 |
| Projection method | Collocated | TVD (Van Leer) | 0.0081 | 0.0034 |
| Projection method | Staggered MAC | TVD (Van Leer) | 0.0070 | 0.0057 |
| Physics-Informed Neural Network | — (mesh-free) | — | 0.1107 | 0.0631 |

All four classical solvers agree with Ghia et al. (1982) to within about 1%, and every TVD variant improves on its upwind counterpart in at least one component — consistent with TVD's reduced numerical diffusion doing what it's supposed to do. The best single result for u is staggered+TVD (0.0070); the best for v is collocated+TVD (0.0034). No single classical combination dominates both components, but all four are close enough to each other that the differences are more a matter of degree than of kind.

The PINN's error is roughly an order of magnitude larger than any classical solver's. This gap did not start there — the first PINN run I trained was considerably worse (max\|u−Ghia\| = 0.2849, max\|v−Ghia\| = 0.2040) — and closing it turned into its own diagnostic exercise, described in the Discussion below.

## Discussion

**The TVD-on-collocated-grid concern from the original note.** I flagged this as an open problem last time: TVD on the collocated grid appeared to interact badly with the pressure Poisson solve. Re-testing it here, the collocated+TVD solver converges cleanly (residual decays monotonically from O(1e-1) to below 1e-6) and achieves the *lowest* v-error of all four classical solvers (0.0034). I don't have a definitive explanation for why this run behaves better than what I described in the original note — it's possible the specific numerical safeguards already in this version of the code (the `eps = 1e-6` regulariser on the TVD slope ratio, which exists specifically to stop the limiter from chattering near-zero gradients close to steady state) are what fixed it, or it's possible the checkerboard modes Rhie & Chow (1983) describe are present but simply don't corrupt the *centreline velocity profile* enough to show up in this error metric even if they'd show up in a direct look at the raw pressure field. I did not go back and inspect the pressure field for checkerboard artefacts in this note, so I'm reporting what the velocity-profile comparison shows and flagging the open question honestly rather than claiming it's fully resolved.

**Why the PINN was so much worse, and what fixed most of it.** The first full training run looked, by its own final loss value, like it had converged — but the printed L-BFGS loss during refinement was going bit-for-bit identical for over 2,700 consecutive optimiser steps, which is not convergence, it's the `strong_wolfe` line search failing to find any improving step. The cause was floating-point precision: the physics-residual loss needs *second*-order derivatives, obtained by differentiating through the network twice (`autograd.grad(create_graph=True)` called twice), and in float32 this accumulates enough round-off noise that L-BFGS's line search runs out of usable precision well short of a real optimum. This exact mechanism — float32 causing L-BFGS to falsely report convergence in PINN training — has been independently documented by Xu et al. (2025); the diagnosis here is a case-study confirmation of a known failure mode applied to a new benchmark, not a new discovery. Switching to float64 training (the original Raissi et al. (2019) reference PINN implementation does this for the same reason) removed the stall completely and cut both errors by roughly half.

A second, compounding issue was loss-term imbalance. The boundary-condition loss weight (`w_b = 100`) was chosen so that boundary constraints would dominate the interior physics residual — standard PINN practice — but by the time training reached a reasonable fit, the *raw*, unweighted boundary-condition error was already about five times smaller than the raw physics-residual error, meaning the 100× weight was making the optimiser spend most of its effort squeezing an already-small boundary error rather than fixing the much larger physics residual that actually encodes the correct nonlinear flow behaviour. This is the loss-imbalance "gradient pathology" documented generally by Wang, Teng & Perdikaris (2021). Reducing the weight in two steps (100 → 20 → 8), based on the measured raw-magnitude ratio rather than a guess, along with fixing an unrelated bug — a hand-written quasi-random sampling function that was both never actually called and, when I tested it directly, turned out to contain a genuine infinite loop for more than about 54 points — brought the final errors down to the PINN row of Table 2.

Even after these fixes, the PINN remains roughly 10× less accurate than any of the classical solvers on this problem. Part of that gap is architectural: at the point these results were recorded, the L-BFGS loss curve had only just started to flatten, suggesting the 21,000-parameter network was approaching the limit of what it could represent, rather than the limit of what the optimiser could reach. A widened network (100 units per hidden layer instead of 64, ~51,000 parameters) has been implemented and locally smoke-tested but not yet trained to completion at the time of writing — see Conclusion.

## Conclusion

Extending the original projection-method solver to all four combinations of {collocated, staggered} × {upwind, TVD} closes the investigation I flagged as ongoing in the original note: TVD convection with a staggered grid and a direct sparse pressure solve works as expected, and — somewhat to my surprise — TVD on the original collocated grid now also converges cleanly and produces the best v-velocity agreement of the four. All four classical solvers agree with Ghia et al. (1982) to within about 1%.

Building a fifth solver on a completely different paradigm — a Physics-Informed Neural Network — let me compare a mesh-free, optimisation-based approach against the classical grid-based ones on identical footing. The PINN's initial results were poor for reasons that had nothing to do with the physics: floating-point precision limits on the optimiser, and a loss-weighting choice that mis-prioritised an already-satisfied constraint over the actual flow physics. Fixing both — training in float64 and rebalancing the loss weights based on measured error magnitudes, plus correcting a broken collocation-sampling routine — cut the PINN's errors substantially, though it still remains behind the classical solvers on this problem at the settings tested. A wider network is queued as the next step; that investigation is, again, ongoing.

## AI Usage Disclosure

Portions of this work — debugging and diagnosing the PINN failure modes described above, implementing the code changes across all five solver notebooks, running and validating the experiments, and drafting this manuscript text — were carried out with the assistance of Claude Code (Anthropic), an AI coding assistant, under the author's direction. All research direction, interpretation, and final editorial decisions are the author's own. Full details are given in the accompanying `AI_USAGE_DISCLOSURE.md`.

## References

Rajaramaduraikarthik, H. R. (2026). *Lid-driven cavity flow: A Python projection method solver validated against Ghia et al. (1982)* [Technical note]. Zenodo. https://doi.org/10.5281/zenodo.20623227

Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. *Journal of Computational Physics, 48*(3), 387–411. https://doi.org/10.1016/0021-9991(82)90058-4

Chorin, A. J. (1968). Numerical solution of the Navier-Stokes equations. *Mathematics of Computation, 22*(104), 745–762. https://doi.org/10.1090/S0025-5718-1968-0242392-2

Patankar, S. V. (1980). *Numerical heat transfer and fluid flow*. Hemisphere Publishing Corporation.

Ferziger, J. H., Perić, M., & Street, R. L. (2020). *Computational methods for fluid dynamics* (4th ed.). Springer. https://doi.org/10.1007/978-3-319-99693-6

Barba, L. A., & Forsyth, G. F. (2018). CFD Python: The 12 steps to Navier-Stokes equations. *Journal of Open Source Education, 2*(16), 21. https://doi.org/10.21105/jose.00021

Harlow, F. H., & Welch, J. E. (1965). Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface. *Physics of Fluids, 8*(12), 2182–2189. https://doi.org/10.1063/1.1761178

van Leer, B. (1974). Towards the ultimate conservative difference scheme II. Monotonicity and conservation combined in a second-order scheme. *Journal of Computational Physics, 14*(4), 361–370. https://doi.org/10.1016/0021-9991(74)90019-9

Sweby, P. K. (1984). High resolution schemes using flux limiters for hyperbolic conservation laws. *SIAM Journal on Numerical Analysis, 21*(5), 995–1011. https://doi.org/10.1137/0721062

Rhie, C. M., & Chow, W. L. (1983). Numerical study of the turbulent flow past an airfoil with trailing edge separation. *AIAA Journal, 21*(11), 1525–1532. https://doi.org/10.2514/3.8284

Raissi, M., Perdikaris, P., & Karniadakis, G. E. (2019). Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations. *Journal of Computational Physics, 378*, 686–707. https://doi.org/10.1016/j.jcp.2018.10.045

Wang, S., Teng, Y., & Perdikaris, P. (2021). Understanding and mitigating gradient flow pathologies in physics-informed neural networks. *SIAM Journal on Scientific Computing, 43*(5), A3055–A3081. https://doi.org/10.1137/20M1318043

Xu, C., Liu, D., Nassereldine, A., & Xiong, J. (2025). FP64 is all you need: Rethinking failure modes in physics-informed neural networks. *Advances in Neural Information Processing Systems (NeurIPS 2025)*. https://arxiv.org/abs/2505.10949
