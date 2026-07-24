from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

from .sampling_rule import SamplingRule, _Rule

@type_check_only
class _SamplingRule(TypedDict):
    version: NotRequired[int]
    default: _Rule
    rules: list[_Rule]

local_sampling_rule: _SamplingRule
SUPPORTED_RULE_VERSION: tuple[int, ...]

class LocalSampler:
    def __init__(self, rules: _SamplingRule = ...) -> None: ...
    def should_trace(self, sampling_req: SamplingRule | None = None) -> bool: ...
    def load_local_rules(self, rules: _SamplingRule) -> None: ...
