import copy
from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.common.parargs import Parargs
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_eval import (
    ResourceEval,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_eval_factory import (
    resource_eval_for,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer import (
    ResourceOutputTransformer,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.eval.resource_output_transformer.transformer_factory import (
    resource_output_transformer_for,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_map.item_reader.reader_config.reader_config_decl import (
    ReaderConfig,
)
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    Resource,
)
from moto.stepfunctions.parser.asl.eval.environment import Environment


class ItemReader(EvalComponent):
    resource_eval: Final[ResourceEval]
    parargs: Final[Optional[Parargs]]
    reader_config: Final[Optional[ReaderConfig]]
    resource_output_transformer: Optional[ResourceOutputTransformer]

    def __init__(
        self,
        resource: Resource,
        parargs: Optional[Parargs],
        reader_config: Optional[ReaderConfig],
    ):
        self.resource_eval = resource_eval_for(resource=resource)
        self.parargs = parargs
        self.reader_config = reader_config

        self.resource_output_transformer = None
        if self.reader_config:
            self.resource_output_transformer = resource_output_transformer_for(
                input_type=self.reader_config.input_type
            )

    @property
    def resource(self):
        return self.resource_eval.resource

    def __str__(self):
        class_dict = copy.deepcopy(self.__dict__)
        del class_dict["resource_eval"]
        class_dict["resource"] = self.resource
        return f"({self.__class__.__name__}| {class_dict})"

    def _eval_body(self, env: Environment) -> None:
        resource_config = None
        if self.reader_config:
            self.reader_config.eval(env=env)
            resource_config = env.stack.pop()

        if self.parargs:
            self.parargs.eval(env=env)
        else:
            env.stack.append({})

        self.resource_eval.eval_resource(env=env)

        if self.reader_config:
            env.stack.append(resource_config)
            self.resource_output_transformer.eval(env=env)
