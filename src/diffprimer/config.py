"""
Configuration loader for DiffPrimer.
"""

from configparser import ConfigParser

from diffprimer.logs import diffprimerLog

logger = diffprimerLog()

DEFAULT_CONFIG = {
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_OPT_SIZE": 22,
    "PRIMER_MIN_SIZE": 20,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_OPT_TM": 64.0,
    "PRIMER_MIN_TM": 60.0,
    "PRIMER_MAX_TM": 68.0,
    "PRIMER_MAX_DIFF_TM": 5.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_OPT_GC": 50.0,
    "PRIMER_MUST_MATCH_FIVE_PRIME": 'SSNNN',
    "PRIMER_GC_CLAMP": 2,
    "PRIMER_MAX_SELF_ANY": 4,
    "PRIMER_MAX_SELF_END": 2,
    "PRIMER_PAIR_MAX_COMPL_ANY": 4,
    "PRIMER_PAIR_MAX_COMPL_END": 2,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_PRODUCT_SIZE_RANGE": [100,300] # default amplicon size range from primer3
}

def load_config(config_file: str | None = None) -> dict:
    """
    Load primer design configuration from an INI file.
    Falls back to default values if no config file is provided.
    Any key defined in the file overwrites the corresponding default.

    Args:
        config_file (str): Path to the configuration INI file. Optional.
    Returns:
        dict: A dictionary with Primer3 configuration parameters.
    """
    global_args = DEFAULT_CONFIG.copy()  # always start from defaults

    if config_file is None:
        logger.info("No configuration file provided. Using default settings.")
        return global_args

    try:
        config = ConfigParser()
        config.read(config_file)
        primer_config = config["PRIMER_SETTINGS"]
        file_args = {key.upper(): eval(value) for key, value in primer_config.items()}
        global_args.update(file_args)  # file values overwrite defaults
        logger.info(f"Configuration loaded from {config_file}")
    except Exception as e:
        raise RuntimeError(logger.error(f"Error loading configuration file: {e}"))

    return global_args
