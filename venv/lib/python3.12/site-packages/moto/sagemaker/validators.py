from typing import Any, Optional


def is_integer_between(
    x: int, mn: Optional[int] = None, mx: Optional[int] = None, optional: bool = False
) -> bool:
    if optional and x is None:
        return True
    try:
        if mn is not None and mx is not None:
            return int(x) >= mn and int(x) < mx
        elif mn is not None:
            return int(x) >= mn
        elif mx is not None:
            return int(x) < mx
        else:
            return True
    except ValueError:
        return False


def is_one_of(x: Any, choices: Any, optional: bool = False) -> bool:
    if optional and x is None:
        return True
    return x in choices
