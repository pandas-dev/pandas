FULL_SAMPLE = 0
ROLLING = 1
EXPANDING = 2

TIME = 0
ENTITY = 1

def _get_cluster_type(cluster_type):
    if cluster_type in (TIME, ENTITY, None):
        return cluster_type

    elif isinstance(cluster_type, basestring):
        cluster_type_up = cluster_type.upper()

        if cluster_type_up == 'ENTITY':
            return ENTITY
        elif cluster_type_up == 'TIME':
            return TIME

    raise Exception('Unrecognized clustering type: %s' % cluster_type)

def _get_window_type(window_type):
    if window_type in (FULL_SAMPLE, ROLLING, EXPANDING):
        return window_type
    elif isinstance(window_type, basestring):
        window_type_up = window_type.upper()

        if window_type_up in ('FULL SAMPLE', 'FULL_SAMPLE'):
            return FULL_SAMPLE
        elif window_type_up == 'ROLLING':
            return ROLLING
        elif window_type_up == 'EXPANDING':
            return EXPANDING

    raise Exception('Unrecognized window type: %s' % window_type)

def banner(text, width=80):
    """

    """
    toFill = width - len(text)

    left = toFill // 2
    right = toFill - left

    return '%s%s%s' % ('-' * left, text, '-' * right)

def f_stat_to_dict(result):
    f_stat, shape, p_value = result

    result = {}
    result['f-stat'] = f_stat
    result['DF X'] = shape[0]
    result['DF Resid'] = shape[1]
    result['p-value'] = p_value

    return result
