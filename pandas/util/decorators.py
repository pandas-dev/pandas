from pandas._tseries import cache_readonly
import warnings

def deprecate(name, alternative):
    alt_name = alternative.func_name
    def wrapper(*args, **kwargs):
        warnings.warn("%s is deprecated. Use %s instead" % (name, alt_name),
                      FutureWarning)
        return alternative(*args, **kwargs)
    return wrapper

