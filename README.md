# Steps to make a fit to LISA mock data

This repository contains source code that can be used to extract thermodynamical parameters characterizing a cosmological phase transition from LISA gravitational-wave mock data by minimizing a fitness function using [evortran](https://gitlab.com/thomas.biekoetter/evortran). Five programs are avaible that are contained in the `app` folder. These programs were used to obtain the results presented in Section 4.2.2 of the evortran paper [arXiv:2507.06082](https://arxiv.org/abs/2507.06082).

1. Clone the repository and navigat to its main folder:
    ```
    git clone https://gitlab.com/thomas.biekoetter/evortran_lisa_fit.git
    cd evortran_lisa_fit
    ```

2. Build the FPM project:
    ```
    source exports/exports_run_gfortran.sh
    fpm build
    ```

3. Run one of the applications:

    ```
    fpm run <application name>
    ```

    The different applications produce the results of the following scenarions described in the evortran paper:

        - `ewpt_general`: Scenario 1: Reconstructing EWPT signal
        - `ewpt_gx_fixed`: Scenario 2: Reconstructing EWPT signal assuming α ≤ 1 and g∗ = 110
        - `ewpt_gx_vw_fixed`: Scenario 3: Reconstructing EWPT signal assuming α ≤ 1, g∗ = 110 and vw = 0.9
        - `ewpt_with_collision_alphalinearprior`: Scenario 4: Reconstructing stronger signal assuming α ≤ 103 (linear prior), g∗ ≤ 150 and vw = 1.
        - `ewpt_with_collision_alphalogprior`: Scenario 4: Reconstructing stronger signal assuming α ≤ 103 (log prior), g∗ ≤ 150 and vw = 1.

    The data is stored in the corresponding subdirectories of the `data` directory. Depending on the number of threads with which the programs are executed, and whether the project is build in debug mode or run mode, the final data can be slightly different compared to the results shown in the plots of the evortran paper.

4. In order to produce the 1 sigma credible region shown in Fig. 12 using `PolyChordLite`, follow the instructions given in `external/PolyChordLite/instructions.md`.

5. In order to produce the 1 sigma credible regions shown in Fig. 12 using `copa`, follow the instructions given in `data/ewpt_gx_vw_fixed_copa/README.md`.

If you use parts of this software in a scientific publication, please cite the following paper:
* [arXiv:2507.06082]: Thomas Biekötter (IFT, Madrid),
    *evortran: a modern Fortran package for genetic algorithms
    with applications from LHC data fitting to
    LISA signal reconstruction*
