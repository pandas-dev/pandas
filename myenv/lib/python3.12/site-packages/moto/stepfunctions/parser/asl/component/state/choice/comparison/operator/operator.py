import abc
from typing import Any

from moto.stepfunctions.parser.asl.eval.environment import Environment
from moto.stepfunctions.parser.utils import SubtypesInstanceManager


class Operator(abc.ABC, SubtypesInstanceManager):
    @staticmethod
    @abc.abstractmethod
    def eval(env: Environment, value: Any) -> None:
        pass
