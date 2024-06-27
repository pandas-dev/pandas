from abc import ABC

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent


class Comparison(EvalComponent, ABC): ...
