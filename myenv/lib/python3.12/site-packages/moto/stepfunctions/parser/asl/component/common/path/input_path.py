import copy
from typing import Final, Optional

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent
from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.asl.utils.json_path import JSONPathUtils


class InputPath(EvalComponent):
    DEFAULT_PATH: Final[str] = "$"

    input_path_src: Final[Optional[str]]

    def __init__(self, input_path_src: Optional[str]):
        self.input_path_src = input_path_src

    def _eval_body(self, env: Environment) -> None:
        if self.input_path_src is None:
            value = dict()
        elif self.input_path_src == InputPath.DEFAULT_PATH:
            value = env.inp
        else:
            value = JSONPathUtils.extract_json(self.input_path_src, env.inp)
        env.stack.append(copy.deepcopy(value))
