"""
Configuration loader for DiffPrimer.
"""

from configparser import ConfigParser
from diffprimer.logs import diffprimerLog

logger = diffprimerLog()

def load_config(config_file: str) -> dict:
    """
    Load primer design configuration from an INI file.
    Args:
        config_file (str): Path to the configuration INI file.
    Returns:
        dict: A dictionary with Primer3 configuration parameters.
    """
    try:
        config = ConfigParser()
        config.read(config_file)
        primer_config = config["PRIMER_SETTINGS"]
        global_args = {key.upper(): eval(value) for key, value in primer_config.items()}
        logger.info(f"Configuration loaded from {config_file}")
    except Exception as e:
        raise RuntimeError(logger.error(f"Error loading configuration file: {e}"))

    return global_args
