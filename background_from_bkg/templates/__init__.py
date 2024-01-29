from .config import config
import yaml


def get_config():
    return yaml.safe_load(config)


__all__ = ["config", "get_config"]
