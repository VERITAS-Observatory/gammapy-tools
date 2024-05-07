from .config import config
import yaml


def get_config() -> dict:
    """Get a template configuration dictionary

    Parameters
    ----------
        None

    Returns
    ----------
        dict                                    - Configuration dictionary


    """
    return yaml.safe_load(config)


__all__ = ["config", "get_config"]
