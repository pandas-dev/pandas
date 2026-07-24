from moto.stepfunctions.parser.api import VariableName, VariableNameList
from moto.stepfunctions.parser.asl.antlr.runtime.ASLIntrinsicParser import (
    ASLIntrinsicParser,
)
from moto.stepfunctions.parser.asl.jsonata.jsonata import (
    VariableReference,
    extract_jsonata_variable_references,
)
from moto.stepfunctions.parser.asl.static_analyser.intrinsic.intrinsic_static_analyser import (
    IntrinsicStaticAnalyser,
)


class VariableNamesIntrinsicStaticAnalyser(IntrinsicStaticAnalyser):
    _variable_names: VariableNameList

    def __init__(self):
        super().__init__()
        self._variable_names = []

    @staticmethod
    def process_and_get(definition: str) -> VariableNameList:
        analyser = VariableNamesIntrinsicStaticAnalyser()
        analyser.analyse(definition=definition)
        return analyser.get_variable_name_list()

    def get_variable_name_list(self) -> VariableNameList:
        return self._variable_names

    def visitFunc_arg_list(self, ctx: ASLIntrinsicParser.Func_arg_listContext) -> None:
        # TODO: the extraction logic is not always in the same order as AWS's
        for child in ctx.children[::-1]:
            self.visit(child)

    def visitFunc_arg_var(self, ctx: ASLIntrinsicParser.Func_arg_varContext) -> None:
        variable_references: set[VariableReference] = (
            extract_jsonata_variable_references(ctx.STRING_VARIABLE().getText())
        )
        for variable_reference in variable_references:
            variable_name: VariableName = variable_reference[1:]
            self._variable_names.append(variable_name)
