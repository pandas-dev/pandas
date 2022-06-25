import importlib


def __getattr__(name):
    if name == "groupby":
        return importlib.import_module("pandas.core.groupby")
    raise AttributeError(f"module 'pandas.core' has no attribute '{name}'")
