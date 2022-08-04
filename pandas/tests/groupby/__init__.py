def get_method_args(name, obj):
    if name in ("nth", "fillna", "take"):
        return (0,)
    if name == "quantile":
        return (0.5,)
    if name == "corrwith":
        return (obj,)
    if name == "tshift":
        return (0, 0)
    return ()
