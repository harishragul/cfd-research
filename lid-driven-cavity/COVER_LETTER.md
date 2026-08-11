Dear Editor,

I am submitting my manuscript, "Effects of Grid Arrangement and Convection-Scheme Discretisation on Projection-Method Accuracy for the Lid-Driven Cavity: A Verification Study," for consideration in the International Journal for Numerical Methods in Fluids.

The manuscript reports a controlled 2×2 verification comparison — collocated versus staggered (MAC) grid arrangement, crossed with first-order upwind versus second-order TVD (Van Leer) convection discretisation — for a projection-method solver applied to the Re = 100 lid-driven cavity, with all four configurations validated against Ghia et al. (1982) using an identical quantitative error metric. This directly addresses your journal's stated interest in verification and validation: while the lid-driven cavity is a long-established benchmark, published studies typically report a single solver configuration, making it difficult to separate how much of the reported accuracy is attributable to grid arrangement versus convection scheme. This manuscript isolates both factors under identical conditions.

Two further contributions are reported: an a priori, quantitative 1D screening study used to select the TVD flux limiter carried into the 2D comparison (rather than a conventional choice), and a re-examination of a collocated-grid TVD convergence difficulty reported in my earlier note, which is found to converge cleanly and to achieve the lowest v-velocity error of the four configurations tested — a result discussed transparently alongside the checkerboard-pressure-mode literature for collocated grids, with the open question this raises explicitly identified rather than overstated as resolved.

The solver code, notebooks, and full numerical output underlying all reported results are publicly available, and are referenced in the manuscript's Data Availability Statement.

I confirm that this manuscript is original, is not under consideration elsewhere, and that all authors have approved its submission. Portions of this work involved AI-assisted coding tools, disclosed in the manuscript's Declaration of AI-Assisted Technologies in accordance with the journal's policy.

Thank you for your consideration.

Sincerely,
Harish Ragul Rajaramaduraikarthik
Indian Institute of Information Technology, Design and Manufacturing, Kancheepuram
harishragulkarthik@gmail.com
