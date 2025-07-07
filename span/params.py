#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyse 1D astronomical spectra.

    1. This software is licensed **for non-commercial use only**.
    2. The source code may be **freely redistributed**, but this license notice must always be included.
    3. Any user who redistributes or uses this software **must properly attribute the original author**.
    4. The source code **may be modified** for non-commercial purposes, but any modifications must be clearly documented.
    5. **Commercial use is strictly prohibited** without prior written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
# Dataclass to store all the GUI parameters to be modified by the user and by the GUI itself

from dataclasses import dataclass, field
import os
import numpy as np
from typing import Union

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_modules import misc

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span.span_modules import misc

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

@dataclass
class SpectraParams:
    # Automatically load result path and create directory structure
    result_path: str = field(default_factory=lambda: SpectraParams.load_result_path())
    result_data: str = field(init=False)  # Initialized in __post_init__

    # Subdirectories
    result_spec_dir: str = field(init=False)
    result_spec: str = field(init=False)
    result_snr_dir: str = field(init=False)
    result_bb_dir: str = field(init=False)
    result_xcorr_dir: str = field(init=False)
    result_vel_disp_dir: str = field(init=False)
    result_ew_data_dir: str = field(init=False)
    result_line_fitting_dir: str = field(init=False)
    result_ppxf_kin_data_dir: str = field(init=False)
    result_ppxf_pop_data_dir: str = field(init=False)
    result_sigma_coeff_dir: str = field(init=False)
    result_plot_dir: str = field(init=False)
    result_long_slit_extract: str = field(init=False)

    spec_names: np.ndarray = field(default_factory=lambda: np.zeros(5))
    spec_names_nopath: np.ndarray = field(default_factory=lambda: np.zeros(5))
    prev_spec: str = ''
    prev_spec_nopath: str = ''
    spectra_list: list = None
    spectra_list_name: list = None
    spectra_number_to_process: int = 0
    spec_names_to_process: list = None
    spec_names_nopath_to_process: list = None
    spec_not_valid: list = None

    # Spectral data
    wavelength: np.ndarray = field(default_factory=lambda: np.array([]))  # Wavelength array
    flux: np.ndarray = field(default_factory=lambda: np.array([]))  # Flux array
    original_wavelength: np.ndarray = field(default_factory=lambda: np.array([]))  # Original wavelength before processing
    original_flux: np.ndarray = field(default_factory=lambda: np.array([]))  # Original flux before processing

    #Continuum model of the Continuum modelling task of the Spectral manipulation panel
    continuum_flux: np.ndarray = field(default_factory=lambda: np.array([]))

    # Processed data
    proc_wavelength: np.ndarray = field(default_factory=lambda: np.array([]))  # Processed wavelength array
    proc_flux: np.ndarray = field(default_factory=lambda: np.array([]))  # Processed flux array

    # Wavelength limits
    wave_limits: np.ndarray = field(init=False)

    #checking parameters
    lambda_units: str = ''
    task_done: int = 0
    task_spec: int = 0
    task_done2: int = 0
    task_spec2: int = 0
    task_analysis: int = 0

   # Task-related parameters of the spectra manipulation panel
    current_order: list = None
    reorder_op: bool = False
    reordered_operations: list = None
    active_operations: list = None


    save_intermediate_spectra: bool = True
    save_final_spectra: bool = False
    not_save_spectra: bool = False

    # Dynamic Cleaning Parameters
    clip_factor: float = 2.5
    sigma_clip_resolution: int = 1600
    sigma_clip_single_vel: bool = True
    sigma_clip_single_value: float = 30.0
    sigma_clip_have_file: bool = False
    sigma_clip_sigma_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "txt_sample_files", "sigma_clip_data.dat"))

    # Denoising Parameters
    moving_average: bool = True
    box_moving_avg: bool = True
    box_moving_avg_size: int = 11
    gauss_moving_avg: bool = False
    gauss_moving_avg_kernel: int = 5
    low_pass_filter: bool = False
    lowpass_cut_off: float = 0.1
    lowpass_order: int = 4
    bandpass_filter: bool = False
    bandpass_lower_cut_off: float = 0.1
    bandpass_upper_cut_off: float = 0.5
    bandpass_order: int = 4

    # Doppler Correction Parameters
    dop_cor_single_shot_vel: float = 0.0
    dop_cor_have_file: bool = False
    dop_cor_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "txt_sample_files", "dopcor_file.dat"))
    dop_cor_single_shot: bool = True
    dop_cor_have_vel: bool = True
    dop_cor_have_z: bool = False

    # Heliocentric Correction Parameters
    helio_have_file: bool = False
    helio_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "txt_sample_files", "file_helio.dat"))
    helio_single_shot: bool = True
    helio_single_shot_location: str = "Paranal"
    helio_single_shot_date: str = "2016-6-4"
    ra_obj: float = 4.88375
    dec_obj: float = 35.0436389

    # Resolution Degradation Parameters
    is_initial_res_r: bool = True
    initial_res_r: int = 1600
    res_degrade_to_r: bool = True
    final_res_r: int = 600
    res_degrade_to_fwhm: bool = False
    final_res_r_to_fwhm: float = 8.4
    is_initial_res_fwhm: bool = False
    initial_res_fwhm: float = 2.51
    final_res_fwhm: float = 8.4
    res_degrade_muse: bool = False
    res_degrade_muse_value: float = 2.51

    # Continuum Subtraction Parameters
    markers_cont_operations: list = field(default_factory=lambda: ["subtract", "divide"])
    cont_math_operation: str = "subtract"  # Initialized in __post_init__
    cont_model_filtering: bool = True
    cont_model_poly: bool = False
    cont_want_to_mask: bool = False
    cont_mask_ranges_str: str = "[(6550, 6660), (4850, 4900),(5850, 5950)]"
    cont_mask_ranges: list[tuple[float, float]] = '' #field(init=False)  # Initialized in __post_init__
    cont_poly_degree: int = 5

    # Blackbody Parameters
    wave1_bb: float = 6000
    wave2_bb: float = 9000
    t_guess: float = 4000.0

    # Cross-Correlation Parameters
    lambda_units_template_crosscorr: str = "a"
    smooth_template_crosscorr: bool = False
    smooth_value_crosscorr: float = 0.0
    low_wave_corr: float = 8400
    high_wave_corr: float = 8800
    is_vel_xcorr: bool = True
    is_z_xcorr: bool = False
    low_vel_corr: float = -1000.0
    high_vel_corr: float = 1000.0
    low_z_corr: float = 0
    high_z_corr: float = 0.1
    lambda_units_template_crosscorr_nm: bool = False
    lambda_units_template_crosscorr_a: bool = True
    lambda_units_template_crosscorr_mu: bool = False

    wave_interval_corr: np.ndarray = field(init=False)
    vel_interval_corr: np.ndarray = field(init=False)
    z_interval_corr: np.ndarray = field(init=False)
    interval_corr: np.ndarray = field(init=False)
    xcorr_limit_wave_range: bool = False
    xcorr_vel_step: float = 5
    xcorr_z_step: float = 0.001

    template_crosscorr: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "templates", "template_emiles.dat"))

    # Sigma Velocity Parameters
    lambda_units_template_sigma: str = "a"
    lambda_units_template_sigma_nm: bool = False
    lambda_units_template_sigma_a: bool = True
    lambda_units_template_sigma_mu: bool = False
    band_cat: bool = True
    band_halpha: bool = False
    band_nad: bool = False
    band_h: bool = False
    band_k: bool = False
    band_custom: bool = False
    resolution_spec: int = 5000
    resolution_template: int = 0
    low_wave_sigma: float = 8400
    high_wave_sigma: float = 8900
    low_wave_cont: float = 8560
    high_wave_cont: float = 8640
    band_sigma: np.ndarray = field(default_factory=lambda: np.array([8440, 8720])) #need to define here because it can be also CaT, Halpha... and then it cannot be linked to any self.parameter
    cont_sigma: np.ndarray = field(default_factory=lambda: np.array([8560, 8640]))

    template_sigma: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "templates", "emiles_template_extended_younger.dat"))

    # Equivalent Width Parameters
    index_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "index_list_sample.txt"))
    have_index_file: bool = False
    single_index: bool = False
    idx_left_blue: float = 8474
    idx_right_blue: float = 8484
    idx_left_red: float = 8563
    idx_right_red: float = 8577
    idx_left_line: float = 8484
    idx_right_line: float = 8513
    index_usr: np.ndarray = field(init=False)

    # LICK/IDS Index Measurements
    lick_ew: bool = True
    lick_constant_fwhm: bool = True
    spec_lick_res_fwhm: float = 3.5
    lick_constant_r: bool = False
    spec_lick_res_r: int = 5000
    radio_lick_sigma_single: bool = False
    radio_lick_sigma_list: bool = False
    radio_lick_sigma_auto: bool = True
    sigma_lick_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "results", "sigma_data_ngc5806_bins_new.dat"))
    sigma_single_lick: int = 100
    correct_ew_sigma: bool = True
    sigma_lick_coeff_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "system_files", "sigma_coeff_lick.dat"))
    lick_index_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "system_files", "lick_indices_angstrom.dat"))
    lick_correct_emission: bool = True
    z_guess_lick_emission: float = 0.0
    lick_ssp_models: list = field(default_factory=lambda: ['Thomas2010', 'xshooter', 'miles', 'smiles'])
    ssp_model: str = 'Thomas2010'  # Impostato in __post_init__
    interp_modes: list = field(default_factory=lambda: ['griddata', 'GPR'])
    interp_model: str = 'GPR'  # Impostato in __post_init__
    stellar_parameters_lick: bool = True
    dop_correction_lick: bool = True

    # Fit Lines Default Values
    emission_line: bool = False
    low_wave_fit: float = 8450
    high_wave_fit: float = 8700
    y0: float = 1.0
    x0: float = 8500
    a: float = -0.8
    sigma: float = 0.5
    m: float = 0.1
    c: float = 1.0
    cat_band_fit: bool = True
    usr_fit_line: bool = False
    wave_interval_fit: np.ndarray = field(init=False)  # in __post_init__
    guess_param: list = field(init=False)  # in __post_init__
    real_cat1: float = 8498
    real_cat2: float = 8542
    real_cat3: float = 8662
    index_ca1: list = field(default_factory=lambda: [8474, 8484, 8563, 8577, 8484, 8513])
    index_ca2: list = field(default_factory=lambda: [8474, 8484, 8563, 8577, 8522, 8562])
    index_ca3: list = field(default_factory=lambda: [8619, 8642, 8700, 8725, 8642, 8682])
    wave_limits_cat: np.ndarray = field(default_factory=lambda: np.array([8440, 8720]))


    # PPXF Kinematics Default Parameters
    kin_stars_templates: list = None
    kin_lam_temp: list = None
    kin_velscale_templates: float = None
    kin_FWHM_gal_cached: Union[float, np.ndarray, None] = None

    wave1_kin: float = 4800
    wave2_kin: float = 5500
    resolution_kin: float = 3.5
    resolution_kin_r: int = 1600
    resolution_kin_muse: bool = False
    sigma_guess_kin: float = 100.0
    redshift_guess_kin: float = 0.0
    constant_resolution_lambda: bool = True
    constant_resolution_r: bool = False
    markers_ppxf_kin: list = field(default_factory=lambda: ['emiles', 'galaxev', 'fsps', 'xshooter'])
    stellar_library_kin: str = 'emiles'
    additive_degree_kin: int = 4
    multiplicative_degree_kin: int = -1
    kin_moments: int = 4
    gas_kin: bool = False
    no_gas_kin: bool = True
    ppxf_kin_noise: float = 0.0163
    kin_best_noise: bool = False
    with_errors_kin: bool = False
    ppxf_kin_preloaded_lib: bool = True
    ppxf_kin_custom_lib: bool = False
    ppxf_kin_lib_folder: str = field(default_factory=lambda: os.path.join(BASE_DIR, "spectralTemplates", "EMILES_BASTI_BASE_KU_FITS"))
    ppxf_kin_custom_temp_suffix: str = '*Eku1.30*.fits'

    ppxf_kin_generic_lib: bool = False
    ppxf_kin_generic_lib_folder: str = field(default_factory=lambda: os.path.join(BASE_DIR, "spectralTemplates", "EMILES_BASTI_BASE_KU_FITS"))
    ppxf_kin_FWHM_tem_generic: float = 2.51
    ppxf_kin_fixed_kin: bool = False

    ppxf_kin_tie_balmer: bool = False
    ppxf_kin_dust_stars: bool = False
    ppxf_kin_dust_gas: bool = False
    ppxf_kin_two_stellar_components: bool = False
    ppxf_kin_age_model1: int = 2
    ppxf_kin_met_model1: int = 0
    ppxf_kin_age_model2: int = 12
    ppxf_kin_met_model2: int = 0
    ppxf_kin_vel_model1: float = 200.0
    ppxf_kin_sigma_model1: float = 200.0
    ppxf_kin_vel_model2: float = 0.0
    ppxf_kin_sigma_model2: float = 50.0
    ppxf_kin_mask_emission: bool = True
    ppxf_kin_have_user_mask: bool = False
    ppxf_kin_mask_ranges_str: str = '[(5180, 5210), (5320, 5380)]'

    ppxf_kin_mc_sim: int = 20
    ppxf_kin_save_spectra: bool = True


    # PPXF Stellar Population Parameters
    pop_with_gas: bool = True
    pop_without_gas: bool = False
    wave1_pop: float = 4800
    wave2_pop: float = 5500
    res_pop: float = 3.5
    z_pop: float = 0.0
    sigma_guess_pop: float = 100.0
    fit_components: str = 'with_gas'
    with_errors: bool = False
    regul_err: float = 0.02
    additive_degree: int = -1
    multiplicative_degree: int = 7
    ppxf_pop_tie_balmer: bool = False
    markers_ppxf: list = field(default_factory=lambda: ['emiles', 'galaxev', 'fsps', 'xshooter', 'sMILES'])

    stellar_library: str = 'emiles'

    ppxf_pop_noise: float = 0.0163
    ppxf_min_age: float = 0.0
    ppxf_max_age: float = 16.0
    ppxf_min_met: float = -2.5
    ppxf_max_met: float = 0.8
    age_range_array: np.ndarray = field(init=False)  # in __post_init__
    met_range_array: np.ndarray = field(init=False)  # in __post_init__
    ppxf_pop_preloaded_lib: bool = True
    ppxf_pop_custom_lib: bool = False
    ppxf_pop_lib_folder: str = field(default_factory=lambda: os.path.join(BASE_DIR, "spectralTemplates", "EMILES_BASTI_BASE_KU_FITS"))
    ppxf_pop_custom_npz: bool = False
    ppxf_pop_npz_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "spectralTemplates", "spectra_emiles_9.0.npz"))
    ppxf_pop_mask: bool = False
    ppxf_custom_temp_suffix: str = '*Eku1.30*.fits'
    ppxf_best_param: bool = False
    ppxf_best_noise_estimate: bool = False
    ppxf_frac_chi: float = 0.30
    ppxf_pop_convolve: bool = False
    ppxf_pop_dust_stars: bool = False
    ppxf_pop_dust_gas: bool = False
    ppxf_pop_want_to_mask: bool = False
    ppxf_pop_mask_ranges_str: str = '[(5180, 5210)]'
    ppxf_pop_mask_ranges_str_default: str = field(init=False)  # in __post_init__
    ppxf_pop_mask_ranges: list = field(init=False)  # in __post_init__
    ppxf_pop_mask_ranges_default: list = field(init=False)  # in __post_init__
    ppxf_pop_error_nsim: int = 20
    ppxf_pop_lg_age: bool = True
    ppxf_pop_lg_met: bool = True
    stellar_parameters_lick_ppxf: bool = False
    lick_ssp_models_ppxf: list = field(default_factory=lambda: ['Thomas2010', 'xshooter', 'miles', 'smiles'])
    ssp_model_ppxf: str = 'Thomas2010'
    interp_modes_ppxf: list = field(default_factory=lambda: ['griddata', 'GPR'])
    interp_model_ppxf: str = 'GPR'

    # Sigma Coeff Parameters
    sigma_coeff: bool = False
    stellar_spectra_coeff_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "sample_templates.txt"))
    lambda_units_coeff_nm: bool = False
    lambda_units_coeff_a: bool = True
    lambda_units_coeff_mu: bool = False
    lambda_units_coeff: str = 'a'
    smooth_stellar_sample: bool = False
    smooth_value_sample: int = 0
    same_idx_ew_task: bool = True
    have_index_file_corr: bool = True
    index_file_corr: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "index_list_sample.txt"))
    single_index_corr: bool = False
    idx_left_blue_sigma: float = 8474
    idx_right_blue_sigma: float = 8484
    idx_left_red_sigma: float = 8563
    idx_right_red_sigma: float = 8577
    idx_left_line_sigma: float = 8484
    idx_right_line_sigma: float = 8513

    # Sigma Correction Default Parameters
    sigma_corr: bool = False
    sigma_vel_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "results", "sigma_data.dat"))
    ew_list_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "results", "ew_data.dat"))
    sigma_coeff_file: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "results", "sigma_coeff.dat"))

    # 2D Long-Slit Extraction Default Parameters
    file_path_spec_extr: str = field(default_factory=lambda: os.path.join(BASE_DIR, "example_files", "spectra/NGC5806_image.fits"))
    trace_y_range_str: str = "(300, 600)"
    poly_degree_str: str = "1"
    extract_y_range_str: str = "(300, 600)"
    snr_threshold_str: str = "20"
    pixel_scale_str: str = "0.252"

    # Cube Extraction Default Parameters
    ifs_run_id: str = 'test'
    ifs_input: str = ''
    ifs_output: str = field(init=False)  #in __post_init__
    ifs_redshift: float = 0.008764
    ifs_lfs_data_default: str = 'None'
    ifs_ow_config: bool = False
    ifs_ow_output: bool = False
    ifs_routine_read: list = field(default_factory=lambda: ['MUSE_WFM', 'MUSE_WFMAOE', 'MUSE_WFMAON', 'MUSE_NFM', 'MUSE_NFMAO', 'CALIFA_V500', 'CALIFA_V1200', 'WEAVE_LIFU', 'JWST_NIRSPEC'])
    ifs_routine_read_default: str = 'MUSE_WFM'
    ifs_origin: str = '14,14'
    ifs_lmin_tot: int = 4800
    ifs_lmax_tot: int = 5500
    ifs_lmin_snr_default: int = 4800
    ifs_lmax_snr_default: int = 5500
    ifs_min_snr_mask: int = 0
    ifs_mask: str = 'none'
    ifs_bin_method: str = 'voronoi'
    ifs_target_snr: int = 50
    ifs_covariance: int = 0
    ifs_prepare_method: str = 'default'
    ifs_preloaded_routine: bool = True
    ifs_user_routine: bool = False
    ifs_user_routine_file: str = ''
    ifs_manual_bin: bool = False
    ifs_voronoi: bool = True
    ifs_existing_bin: bool = False
    ifs_existing_bin_folder: str = ''

    # Spectra Pre-processing Default Parameters
    cropping_spectrum: bool = False
    cropping_low_wave: float = 4800
    cropping_high_wave: float = 5500
    sigma_clipping: bool = False
    wavelet_cleaning: bool = False
    sigma_wavelets: float = 0.02
    wavelets_layers: int = 3
    filter_denoise: bool = False
    dop_cor: bool = False
    helio_corr: bool = False

    # Spectra Processing Frame Default Parameters
    rebinning: bool = False
    rebinning_log: bool = False
    rebinning_linear: bool = True
    rebin_step_pix: float = 0.02
    rebin_step_sigma: int = 60
    degrade: bool = False
    normalize_wave: bool = False
    norm_wave: float = 5000
    sigma_broad: bool = False
    sigma_to_add: float = 0.0
    add_noise: bool = False
    noise_to_add: float = 10.0
    continuum_sub: bool = False

    # Math Frame Default Parameters
    average_all: bool = False
    norm_and_average: bool = False
    do_nothing: bool = True
    sum_all: bool = False
    normalize_and_sum_all: bool = False
    use_for_spec_an: bool = False
    subtract_normalized_avg: bool = False
    subtract_normalized_spec: bool = False
    spectra_to_subtract: str = 'Spectrum to subtract'
    add_pedestal: bool = False
    pedestal_to_add: float = 0.0
    multiply: bool = False
    multiply_factor: float = 1.0
    derivatives: bool = False

    # Spectra Manipulation - Task Management
    active_operations: list = field(default_factory=list)

    # Variables to prevent crashes in case some loaded spectra are not valid
    spectra_number: int = 0
    fatal_condition: int = 0


    def __post_init__(self):
        """Initialize dependent attributes after object creation."""
        self.result_data = misc.create_result_structure(self.result_path)

        # Define subdirectories
        self.result_spec_dir = os.path.join(self.result_data, 'processed_spectra')
        self.result_spec = os.path.join(self.result_spec_dir, '')
        self.result_snr_dir = os.path.join(self.result_data, 'SNR')
        self.result_bb_dir = os.path.join(self.result_data, 'planck_black_body_fitting')
        self.result_xcorr_dir = os.path.join(self.result_data, 'cross-correlation')
        self.result_vel_disp_dir = os.path.join(self.result_data, 'velocity_dispersion')
        self.result_ew_data_dir = os.path.join(self.result_data, 'line-strength_analysis')
        self.result_line_fitting_dir = os.path.join(self.result_data, 'line_fitting')
        self.result_ppxf_kin_data_dir = os.path.join(self.result_data, 'stars_and_gas_kinematics')
        self.result_ppxf_pop_data_dir = os.path.join(self.result_data, 'stellar_populations_and_sfh')
        self.result_sigma_coeff_dir = os.path.join(self.result_data, 'line-strength_sigma_coefficients')
        self.result_plot_dir = os.path.join(self.result_data, 'plots')
        self.result_long_slit_extract = os.path.join(self.result_data, 'longslit_extracted')

        # Initialize computed fields

        try:
            self.wave_limits = np.array([self.wavelength[0], self.wavelength[-1]])
            self.wave_interval_corr = np.array([self.low_wave_corr, self.high_wave_corr])
            self.vel_interval_corr = np.array([self.low_vel_corr, self.high_vel_corr])
            self.z_interval_corr = np.array([self.low_z_corr, self.high_z_corr])
            self.interval_corr = np.array([self.low_z_corr, self.high_z_corr])

            self.cont_mask_ranges = eval(self.cont_mask_ranges_str)
            self.index_usr = np.array([self.idx_left_blue, self.idx_right_blue, self.idx_left_red, self.idx_right_red, self.idx_left_line, self.idx_right_line]).T

            self.wave_interval_fit = np.array([self.low_wave_fit, self.high_wave_fit])
            self.guess_param = [self.y0, self.x0, self.a, self.sigma, self.m, self.c]
            self.age_range_array = np.array([self.ppxf_min_age, self.ppxf_max_age])
            self.met_range_array = np.array([self.ppxf_min_met, self.ppxf_max_met])

        except Exception:
            self.wave_limits = np.array([])
            self.wave_interval_corr = np.array([])
            self.vel_interval_corr = np.array([])
            self.z_interval_corr = np.array([])
            self.interval_corr = np.array([])
            self.cont_mask_ranges = self.cont_mask_ranges_str
            self.index_usr = np.array([])

            self.wave_interval_fit = np.array([])
            self.guess_param = []

        # PPXF kin Masking
        try:
            self.ppxf_kin_mask_ranges = eval(self.ppxf_kin_mask_ranges_str)
        except Exception:
            self.ppxf_kin_mask_ranges = self.ppxf_kin_mask_ranges_str

        # PPXF Stellar Population Masking
        try:
            self.ppxf_pop_mask_ranges_str_default = self.ppxf_pop_mask_ranges_str
            self.ppxf_pop_mask_ranges = eval(self.ppxf_pop_mask_ranges_str)
            self.ppxf_pop_mask_ranges_default = self.ppxf_pop_mask_ranges
        except Exception:
            self.ppxf_pop_mask_ranges_str_default = self.ppxf_pop_mask_ranges_str
            self.ppxf_pop_mask_ranges = self.ppxf_pop_mask_ranges_str
            self.ppxf_pop_mask_ranges_default = self.ppxf_pop_mask_ranges

        # Age and metallicity range arrays
        self.age_range_array = np.array([self.ppxf_min_age, self.ppxf_max_age])
        self.met_range_array = np.array([self.ppxf_min_met, self.ppxf_max_met])

        # Cube extraction output path
        self.ifs_output = self.result_data + '/'


    @staticmethod
    def load_result_path():
        """Load result_path from config.json, or ask the user if missing."""
        config_file = os.path.join(BASE_DIR, "system_files", "config.json")
        config_folder = misc.load_config(config_file)

        if "result_path" not in config_folder or not os.path.exists(config_folder["result_path"]):
            result_path = misc.ask_user_for_result_path()
            if result_path:
                config_folder["result_path"] = result_path
                misc.save_config_folder(config_folder, config_file)
            else:
                sg.popup("No path selected, the program will close")
                exit()

        return config_folder["result_path"]
