from .enums import *

default_fear_function = FearFunctions.FearDisabled.value

default_initial_conditions = ([{
    PERSON_INDEX: 0,
    CONTRACTION_TIME: 0,
    INFECTION_STATUS: InfectionStatus.Contraction
}])

default_stop_simulation_threshold = 10000

default_detection_status = DetectionStatus.NotDetected.value

default_quarantine_status = QuarantineStatus.NoQuarantine.value

default_detection_mild_proba = 0.5

default_distribution = {
    DISTRIBUTION: 'poisson'
}

default_case_severity_distribution = {
    ASYMPTOMATIC: 0.006,
    MILD: 0.809,
    SEVERE: 0.138,
    CRITICAL: 0.047
}

default_disease_times_distributions = {
    T0: default_distribution,
    T1: default_distribution,
    T2: default_distribution
}

default_disease_progression = {
    NOT_DETECTED: default_disease_times_distributions,
    DETECTED: default_disease_times_distributions
}

default_fear_factor = {
    FEAR_FUNCTION: default_fear_function,
    SCALE_FACTOR: 10,
    LOC_FACTOR: 0,
    LIMIT_VALUE: 0.5,
    DETECTED_MULTIPLIER: 0,
    DEATHS_MULTIPLIER: 1
}

default_fear_factors = {DEFAULT: default_fear_factor}

default_experiment_id = 'experiment'

default_transmission_probabilities = {
        transmission_way: 1.0 for transmission_way in KernelType.map()
}

default_import_intensity = {
    FUNCTION: NO_IMPORT
}

default_start_time = 0.0

default_output_root_dir = 'outputs'

default_save_input_data = True

default_death_probability = {
    ASYMPTOMATIC: 0.0,
    MILD: 0.0,
    SEVERE: 0.0,
    CRITICAL: 0.49
}

default_random_seed = 42

default_log_outputs = True
default_max_time = 1000 #float('inf')

default_icu_availability = 100
default_hospital_beds_availability = 5000
default_med_personnel_availability = 400

default_log_time_freq = 1.0

default_serial_interval = {
    MIN_TIME: 0.0,
    MAX_TIME: 30.0
}

default_turn_on_detection = True

default_average_infectivity_time_constant_kernel = 2.339

default_save_expected_severity = False

default_move_zero_time_according_to_detected = False
default_number_of_detected_at_zero_time = 0

default_use_today_mark = False
default_today_offset = 0
default_laid_curve = {}

default_plot_xlim_right = None
default_plot_ylim_top = None
default_plot_xlim_left = 0
default_plot_ylim_bottom = 0

default_plot_xlim_cut_right = None
default_plot_ylim_cut_top = None
default_plot_xlim_cut_left = None
default_plot_ylim_cut_bottom = None

default_enable_visualization = False

default_r_out_schedule = []

default_enable_additional_logs = False
default_reuse_expected_case_severities = False
default_reuse_time_distribution_realizations = False

default_stop_simulation_threshold_type = PREVALENCE

default_old_implementation_for_household_kernel = False

default_constant_age_setup = None

default_inter_age_contacts = False

defaults = {
    INITIAL_CONDITIONS: default_initial_conditions,
    STOP_SIMULATION_THRESHOLD: default_stop_simulation_threshold,
    DISEASE_PROGRESSION: default_disease_progression,
    CASE_SEVERITY_DISTRIBUTION: default_case_severity_distribution,
    OUTPUT_ROOT_DIR: default_output_root_dir,
    EXPERIMENT_ID: default_experiment_id,
    SAVE_INPUT_DATA: default_save_input_data,
    TRANSMISSION_PROBABILITIES: default_transmission_probabilities,
    FEAR_FACTORS: default_fear_factors,
    IMPORT_INTENSITY: default_import_intensity,
    START_TIME: default_start_time,
    DEATH_PROBABILITY: default_death_probability,
    RANDOM_SEED: default_random_seed,
    MAX_TIME: default_max_time,
    ICU_AVAILABILITY: default_icu_availability,
    HOSPITAL_BEDS_AVAILABILITY: default_hospital_beds_availability,
    MED_PERSONNEL_AVAILABILITY: default_med_personnel_availability,
    LOG_TIME_FREQ: default_log_time_freq,
    LOG_OUTPUTS: default_log_outputs,
    SERIAL_INTERVAL: default_serial_interval,
    DETECTION_MILD_PROBA: default_detection_mild_proba,
    TURN_ON_DETECTION: default_turn_on_detection,
    AVERAGE_INFECTIVITY_TIME_CONSTANT_KERNEL: default_average_infectivity_time_constant_kernel,
    SAVE_EXPECTED_SEVERITY: default_save_expected_severity,
    MOVE_ZERO_TIME_ACCORDING_TO_DETECTED: default_move_zero_time_according_to_detected,
    NUMBER_OF_DETECTED_AT_ZERO_TIME: default_number_of_detected_at_zero_time,
    USE_TODAY_MARK: default_use_today_mark,
    TODAY_OFFSET: default_today_offset,
    LAID_CURVE: default_laid_curve,
    PLOT_XLIM_CUT_LEFT: default_plot_xlim_cut_left,
    PLOT_XLIM_CUT_RIGHT: default_plot_xlim_cut_right,
    PLOT_XLIM_LEFT: default_plot_xlim_left,
    PLOT_XLIM_RIGHT: default_plot_xlim_right,
    PLOT_YLIM_CUT_BOTTOM: default_plot_ylim_cut_bottom,
    PLOT_YLIM_CUT_TOP: default_plot_ylim_cut_top,
    PLOT_YLIM_BOTTOM: default_plot_ylim_bottom,
    PLOT_YLIM_TOP: default_plot_ylim_top,
    STOP_SIMULATION_THRESHOLD_TYPE: default_stop_simulation_threshold_type,
    ENABLE_VISUALIZATION: default_enable_visualization,
    R_OUT_SCHEDULE: default_r_out_schedule,
    ENABLE_ADDITIONAL_LOGS: default_enable_additional_logs,
    REUSE_EXPECTED_CASE_SEVERITIES: default_reuse_expected_case_severities,
    REUSE_TIME_DISTRIBUTION_REALIZATIONS: default_reuse_time_distribution_realizations,
    OLD_IMPLEMENTATION_FOR_HOUSEHOLD_KERNEL: default_old_implementation_for_household_kernel,
    CONSTANT_AGE_SETUP: default_constant_age_setup,
}

default_age_induced_fatality_rates = [(0, 20, 0.002), (20, 40, 0.002), (40, 50, 0.004), (50, 60, 0.013),
                                      (60, 70, 0.036), (70, 80, 0.08), (80, 200, 0.148)]

default_age_cohorts_with_descriptions = [(0, 20, '[0-19]'), (20, 40, '[20-39]'), (40, 50, '[40-49]'),
                                         (50, 60, '[50-59]'), (60, 70, '[60-69]'), (70, 80, '[70-79]'),
                                         (80, 200, '[80+]')]