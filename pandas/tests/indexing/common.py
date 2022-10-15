""" common utilities """


def _mklbl(prefix, n):
    return [f"{prefix}{i}" for i in range(n)]


def check_result(obj, method, key, axes=None, fails=None):
    if axes is None:
        axes = [0, 1]
    else:
        assert axes in [0, 1]
        axes = [axes]

    for ax in axes:
        if ax < obj.ndim:
            # create a tuple accessor
            axes = [slice(None)] * obj.ndim
            axes[ax] = key
            axified = tuple(axes)
            try:
                getattr(obj, method).__getitem__(axified)
            except (IndexError, TypeError, KeyError) as detail:
                # if we are in fails, the ok, otherwise raise it
                if fails is not None:
                    if isinstance(detail, fails):
                        return
                raise
