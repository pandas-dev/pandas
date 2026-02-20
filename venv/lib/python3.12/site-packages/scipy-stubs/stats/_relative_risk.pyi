from dataclasses import dataclass

from ._common import ConfidenceInterval

@dataclass
class RelativeRiskResult:
    relative_risk: float
    exposed_cases: int
    exposed_total: int
    control_cases: int
    control_total: int
    def confidence_interval(self, /, confidence_level: float = 0.95) -> ConfidenceInterval: ...

def relative_risk(exposed_cases: int, exposed_total: int, control_cases: int, control_total: int) -> RelativeRiskResult: ...
