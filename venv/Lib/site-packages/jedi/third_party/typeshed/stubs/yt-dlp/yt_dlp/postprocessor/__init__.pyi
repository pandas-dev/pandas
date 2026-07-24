from collections.abc import Mapping

from .common import PostProcessor

_default_pps: Mapping[str, type[PostProcessor]]

def get_postprocessor(key: str) -> type[PostProcessor]: ...
