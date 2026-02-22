"""utilities for working with optional depednencies."""

import json
import os
import warnings
from functools import lru_cache
from typing import Optional


@lru_cache(100)
def has_optional_dependency(importable: str, hint: Optional[str] = None) -> bool:
    """whether a given optional dependency is even installed, with an optional hint"""
    env_var = f"""JUPYTERLITE_NO_{importable.upper().replace('.','_')}"""

    if env_var in os.environ and bool(json.loads(os.environ[env_var])):
        return False

    try:
        __import__(importable)
        return True
    except Exception as error:
        if hint:
            warnings.warn(hint.format(error=error), stacklevel=2)
        return False
