from typing import Set

from moto.stepfunctions.parser.asl.antlr.runtime.ASLParser import ASLParser
from moto.stepfunctions.parser.asl.component.state.exec.state_task.service.resource import (
    ActivityResource,
    Resource,
    ServiceResource,
)
from moto.stepfunctions.parser.asl.component.state.state_type import StateType
from moto.stepfunctions.parser.asl.parse.test_state.asl_parser import (
    TestStateAmazonStateLanguageParser,
)
from moto.stepfunctions.parser.asl.static_analyser.static_analyser import StaticAnalyser


class TestStateStaticAnalyser(StaticAnalyser):
    _SUPPORTED_STATE_TYPES: Set[StateType] = {
        StateType.Task,
        StateType.Pass,
        StateType.Wait,
        StateType.Choice,
        StateType.Succeed,
        StateType.Fail,
    }

    def analyse(self, definition) -> None:
        _, parser_rule_context = TestStateAmazonStateLanguageParser.parse(definition)
        self.visit(parser_rule_context)

    def visitState_type(self, ctx: ASLParser.State_typeContext) -> None:
        state_type_value: int = ctx.children[0].symbol.type
        state_type = StateType(state_type_value)
        if state_type not in self._SUPPORTED_STATE_TYPES:
            raise ValueError(
                f"Unsupported state type for TestState runs '{state_type}'."
            )

    def visitResource_decl(self, ctx: ASLParser.Resource_declContext) -> None:
        resource_str: str = ctx.keyword_or_string().getText()[1:-1]
        resource = Resource.from_resource_arn(resource_str)

        if isinstance(resource, ActivityResource):
            raise ValueError(
                f"ActivityResources are not supported for TestState runs {resource_str}."
            )

        if isinstance(resource, ServiceResource):
            if resource.condition is not None:
                raise ValueError(
                    f"Service integration patterns are not supported for TestState runs {resource_str}."
                )
