from typing import Any, Dict, List

from moto.stepfunctions.parser.asl.component.common.assign.assign_decl_binding import (
    AssignDeclBinding,
)
from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment


class AssignDecl(EvalComponent):
    declaration_bindings: List[AssignDeclBinding]

    def __init__(self, declaration_bindings: List[AssignDeclBinding]):
        super().__init__()
        self.declaration_bindings = declaration_bindings

    def _eval_body(self, env: Environment) -> None:
        declarations: Dict[str, Any] = dict()
        for declaration_binding in self.declaration_bindings:
            declaration_binding.eval(env=env)
            binding: Dict[str, Any] = env.stack.pop()
            declarations.update(binding)
        for identifier, value in declarations.items():
            env.variable_store.set(variable_identifier=identifier, variable_value=value)
