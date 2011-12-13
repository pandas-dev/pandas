def _get_cluster_type(cluster_type):
    cluster_type = _WINDOW_TYPES.get(cluster_type, cluster_type)
    if cluster_type is None:
        return cluster_type

    cluster_type_up = cluster_type.upper()

    if cluster_type_up == 'ENTITY':
        return 'entity'
    elif cluster_type_up == 'TIME':
        return 'time'
    else:  # pragma: no cover
        raise Exception('Unrecognized cluster type: %s' % cluster_type)

_CLUSTER_TYPES = {
    0 : 'time',
    1 : 'entity'
}

_WINDOW_TYPES = {
    0 : 'full_sample',
    1 : 'rolling',
    2 : 'expanding'
}


def _get_window_type(window_type):
    window_type = _WINDOW_TYPES.get(window_type, window_type)
    window_type_up = window_type.upper()

    if window_type_up in ('FULL SAMPLE', 'FULL_SAMPLE'):
        return 'full_sample'
    elif window_type_up == 'ROLLING':
        return 'rolling'
    elif window_type_up == 'EXPANDING':
        return 'expanding'
    else:  # pragma: no cover
        raise Exception('Unrecognized window type: %s' % window_type)

def banner(text, width=80):
    """

    """
    toFill = width - len(text)

    left = toFill // 2
    right = toFill - left

    return '%s%s%s' % ('-' * left, text, '-' * right)
