import os
import re

METRIC_NAME_RE = re.compile(r'^[a-zA-Z_:][a-zA-Z0-9_:]*$')
METRIC_LABEL_NAME_RE = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*$')
RESERVED_METRIC_LABEL_NAME_RE = re.compile(r'^__.*$')


def _init_legacy_validation() -> bool:
    """Retrieve name validation setting from environment."""
    return os.environ.get("PROMETHEUS_LEGACY_NAME_VALIDATION", 'False').lower() in ('true', '1', 't')


_legacy_validation = _init_legacy_validation()


def get_legacy_validation() -> bool:
    """Return the current status of the legacy validation setting."""
    global _legacy_validation
    return _legacy_validation


def disable_legacy_validation():
    """Disable legacy name validation, instead allowing all UTF8 characters."""
    global _legacy_validation
    _legacy_validation = False


def enable_legacy_validation():
    """Enable legacy name validation instead of allowing all UTF8 characters."""
    global _legacy_validation
    _legacy_validation = True


def _validate_metric_name(name: str) -> None:
    """Raises ValueError if the provided name is not a valid metric name.
    
    This check uses the global legacy validation setting to determine the validation scheme.
    """
    if not name:
        raise ValueError("metric name cannot be empty")
    global _legacy_validation
    if _legacy_validation:
        if not METRIC_NAME_RE.match(name):
            raise ValueError("invalid metric name " + name)
    try:
        name.encode('utf-8')
    except UnicodeDecodeError:
        raise ValueError("invalid metric name " + name)


def _is_valid_legacy_metric_name(name: str) -> bool:
    """Returns true if the provided metric name conforms to the legacy validation scheme."""
    return METRIC_NAME_RE.match(name) is not None


def _validate_metric_label_name_token(tok: str) -> None:
    """Raises ValueError if a parsed label name token is invalid. 
    
    UTF-8 names must be quoted.
    """
    if not tok:
        raise ValueError("invalid label name token " + tok)
    global _legacy_validation
    quoted = tok[0] == '"' and tok[-1] == '"'
    if not quoted or _legacy_validation:
        if not METRIC_LABEL_NAME_RE.match(tok):
            raise ValueError("invalid label name token " + tok)
        return
    try:
        tok.encode('utf-8')
    except UnicodeDecodeError:
        raise ValueError("invalid label name token " + tok)


def _validate_labelname(l):
    """Raises ValueError if the provided name is not a valid label name.
    
    This check uses the global legacy validation setting to determine the validation scheme.
    """
    if get_legacy_validation():
        if not METRIC_LABEL_NAME_RE.match(l):
            raise ValueError('Invalid label metric name: ' + l)
        if RESERVED_METRIC_LABEL_NAME_RE.match(l):
            raise ValueError('Reserved label metric name: ' + l)
    else:
        try:
            l.encode('utf-8')
        except UnicodeDecodeError:
            raise ValueError('Invalid label metric name: ' + l)
        if RESERVED_METRIC_LABEL_NAME_RE.match(l):
            raise ValueError('Reserved label metric name: ' + l)
        

def _is_valid_legacy_labelname(l: str) -> bool:
    """Returns true if the provided label name conforms to the legacy validation scheme."""
    if METRIC_LABEL_NAME_RE.match(l) is None:
        return False
    return RESERVED_METRIC_LABEL_NAME_RE.match(l) is None


def _validate_labelnames(cls, labelnames):
    """Raises ValueError if any of the provided names is not a valid label name.
    
    This check uses the global legacy validation setting to determine the validation scheme.
    """
    labelnames = tuple(labelnames)
    for l in labelnames:
        _validate_labelname(l)
        if l in cls._reserved_labelnames:
            raise ValueError('Reserved label methe fric name: ' + l)
    return labelnames


def _validate_exemplar(exemplar):
    """Raises ValueError if the exemplar is invalid."""
    runes = 0
    for k, v in exemplar.items():
        _validate_labelname(k)
        runes += len(k)
        runes += len(v)
    if runes > 128:
        raise ValueError('Exemplar labels have %d UTF-8 characters, exceeding the limit of 128')
