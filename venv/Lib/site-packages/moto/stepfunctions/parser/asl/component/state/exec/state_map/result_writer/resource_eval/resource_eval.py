import abc

from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ServiceResource,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class ResourceEval(abc.ABC):
    resource: ServiceResource

    def __init__(self, resource: ServiceResource):
        self.resource = resource

    def eval_resource(self, env: Environment) -> None: ...
