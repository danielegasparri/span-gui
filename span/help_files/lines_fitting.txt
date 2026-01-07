SPAN: SPectral ANalysis software V7.4
Daniele Gasparri, January 2026

# Line(s) Fitting

Line(s) fitting in SPAN measures the properties of spectral lines (centers, widths, and integrated fluxes), plus uncertainties, within a user-defined wavelength window or for the Ca II Triplet (CaT). Typical uses include: gas emission-line fluxes and kinematics proxies, stellar absorption features (e.g., CaT strength/width), and quick checks of blended features in binned spectra.

SPAN offers two modes:

- CaT fitting (automatic): A specialized, fully automatic fit for the Ca II triplet around 8500-8700 A, using three gaussian functions and a line with slope for the continuum. In this case, your spectra must be rest-frame corrected and/or Doppler corrected.
- Generic lines fitting (custom): A flexible fit for one or more lines in any window you choose, with robust continuum handling and optional bootstrap errors.

Both modes print a compact summary in the terminal window and can optionally save a diagnostic plot per spectrum. 

## CaT mode (automatic)
This option selects a fixed CaT window (default 8440-8720 A), builds a robust local continuum (linear baseline) and normalizes the spectrum. Therefore, it fits three Gaussians (one per triplet line) with sensible constraints and automatic component ordering to the nominal air wavelengths (8498.02, 8542.09, 8662.14 A). For every CaT line, it estimates the integrated flux in physical units (scaled back from the normalized model), computes the equivalent width (EW) from the model (flux/baseline at the line center), and line widths (sigma in A and km/s). Finally, it runs a bootstrap (N=100) by default to return 1sigma errors on centers, widths, fluxes and EWs.

### Outputs

For each spectrum, you get: central wavelength (A), sigma (A), sigma (km/s), flux, EW, and their errors, for each triplet line. A compact per-spectrum summary is printed to the terminal window. In 'Process all' mode, and if 'Save spectra analysis plots' is enabled, a plot with spectrum, model and residuals is saved.

The 'Process all' mode writes a CaT results file (space-separated .dat) with one row per line per spectrum (columns described below).


## Custom line(s) fitting

This mode is designed for any window and any reasonable line shape (Gaussian or Lorentzian), with automatic peak detection and multi-component fitting.

**Parameters you can set:**  

- Window: lambda min, lambda max: The wavelength interval (in A) to consider for line(s) detection and continuum level estimation.
- Line type: auto | emission | absorption. Which kind of line(s) you are trying to identify and fit. The option 'auto' works fine, but for better results, please specify if the lines are in emission or absorption. 
- Peak threshold (sigma): Minimum prominence (in units of the local noise estimate) for a peak to be considered. **Tip:** start with 10-15 for strong emission lines, 3-4 for absorption lines and high (50-100) S/N spectra.
- Components: auto | fixed (+ N). Auto: detect peaks and fit up to Max components. Fixed: fit exactly N components (1-3) placed from detected/seed peaks.
- Max components: If 'Components' is set to 'auto', this will tell SPAN the maximum number of lines to consider from the peak detection results, starting from the strongest. You can fit at maximum 5 lines. **Tip:** keep small (2-3) in narrow windows to avoid over-fitting.
- Profile: gauss | lorentz. Choose the analytic form for components. **Tip:** Gaussians are standard for most nebular and stellar features; Lorentzians can help with extended wings (be cautious with interpretation).
- Resolution (FWHM): The instrumental resolution of your spectra. This parameter is optional. If provided, SPAN will report intrinsic widths for Gaussians via sigma_intr^2 = sigma_obs^2 - sigma_inst^2 (truncated at zero). If omitted, only observed sigma are reported.
- Continuum modelling: Auto (recommended): robust, slope-aware baseline selection per sign. Flat: constant continuum (use for very small windows). Linear: straight line (good default if continuum shows a mild slope). Binned percentile: rolling/large-scale background via percentile in bins (set bin width A and percentiles). Poly2: quadratic baseline (use sparingly; risk of absorbing broad lines).
- Additional settings:
    - Percentile (emission) (default 15): percentile used when expecting emission (lower is safer).
    - Percentile (absorption) (default 85): percentile used for absorption (higher is safer).
    - Bin width (default 50): bin width for the binned baseline.


**Uncertainties**  

Uncertainties estimation is optional and is carried with bootstrap simulations: enable and set N boot (e.g., 50-200). Returns errors on line centers, widths, and fluxes. If bootstrap is off, asymptotic (covariance-based) errors are returned when reliable.

### Outputs

For each fitted component, SPAN returns: Center (A) with error, width sigma (A) and sigma (km/s) with errors, Flux (integrated, physical units) with error, derived from the normalized model and local baseline scale.
