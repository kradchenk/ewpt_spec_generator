This code generates GW spectra coming from electroweak phase transitions according to Lisa Cosmology WG recommendations given in [2403.03723]. The inputs are the macroscopic parameters: Tn, alpha, beta/H, g, vw

Installation:
    ```
    https://github.com/kradchenk/ewpt_spec_generator.git
    cd ewpt_spec_generator
    source exports/exports_run_gfortran.sh
    fpm build
    ```

As an example one can run:

    ```
    fpm run ewpt_signal
    ```

The existing applications are:

        - `ewpt_signal`: Generates the EWPT spectrum corresponding to some parameters, by default includes only sound wave contribution and turbulence (eps = 0.5)
        - `generate_samples`:  Generates a number of random signals within some specified bounds for the macro parameters, it stores the signal + LISA noise generated spectrum. The binning is logarithmic.
        - `generate_samples_var`: Same as before but stores also the variance per bin.
        - `generate_samples_bin`: Same as before but using the binning algorithm decribed in the first paragraph of page 15 of Ref 1906.09244.


All data is stored `data` directory, the main scv contains the parameters of each curve and the generated spectrum, another csv contains the frequencies.
Please check the example in the python directory to treat the data.