import importlib


_delayed_symbols = {
    "Dict": ".typeddict",
    "List": ".typedlist",
}


def __getattr__(name):
    # Uses PEP-562 but requires python>3.6
    if name in _delayed_symbols:
        modpath = _delayed_symbols[name]
        mod = importlib.import_module(modpath, __name__)
        return getattr(mod, name)
    else:
        try:
            return importlib.import_module(f".{name}", __name__)
        except ModuleNotFoundError:
            raise AttributeError
