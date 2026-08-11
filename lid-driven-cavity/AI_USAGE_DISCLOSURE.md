# AI Usage Disclosure

This work used Claude Code (Anthropic), an AI coding assistant, in the following ways:

- **Debugging and root-cause diagnosis** of the PINN solver's poor initial agreement with the Ghia et al. (1982) benchmark, including identifying that L-BFGS was silently stalling under float32 precision, that a boundary-condition loss weight was overwhelming the physics-residual gradient, and that a hand-written Halton low-discrepancy sampling function was unused dead code containing an unrelated infinite-loop bug.
- **Implementing code fixes and improvements** across all five solver notebooks: float64 training precision, loss-weight rebalancing, corrected quasi-random collocation sampling with boundary-layer clustering, network capacity increases, and a standardized quantitative error metric (interpolation against the Ghia tabulated points) added consistently across all five methods for direct comparison.
- **Running and validating experiments**, including local smoke tests of each change before handing off to GPU execution, and executing the four classical finite-difference solvers to steady state to obtain the results reported here.
- **Drafting manuscript text** for this document, based on the author's original preprint (Zenodo, DOI 10.5281/zenodo.20623227) as a style and structure reference, and on the real, author-supervised experimental results described above.

All research direction, interpretation of results, and final editorial decisions were made by the author. The author reviewed and take full responsibility for the code, the results, and the manuscript content.
