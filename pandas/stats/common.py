
_WINDOW_TYPES = {
    0: 'full_sample',
    1: 'rolling',
    2: 'expanding'
}
# also allow 'rolling' as key
_WINDOW_TYPES.update((v, v) for k,v in _WINDOW_TYPES.items())
_ADDITIONAL_CLUSTER_TYPES = set(("entity", "time"))

def _get_cluster_type(cluster_type):
    # this was previous behavior
    if cluster_type is None:
        return cluster_type
    try:
        return _get_window_type(cluster_type)
    except ValueError:
        final_type = str(cluster_type).lower().replace("_", " ")
        if final_type in _ADDITIONAL_CLUSTER_TYPES:
            return final_type
        raise ValueError('Unrecognized cluster type: %s' % cluster_type)

def _get_window_type(window_type):
    # e.g., 0, 1, 2
    final_type = _WINDOW_TYPES.get(window_type)
    # e.g., 'full_sample'
    final_type = final_type or _WINDOW_TYPES.get(str(window_type).lower().replace(" ", "_"))
    if final_type is None:
        raise ValueError('Unrecognized window type: %s' % window_type)
    return final_type

def banner(text, width=80):
    """

    """
    toFill = width - len(text)

    left = toFill // 2
    right = toFill - left

    return '%s%s%s' % ('-' * left, text, '-' * right)
