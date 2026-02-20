import abc

from moto.stepfunctions.parser.asl.component.eval_component import EvalComponent


class ResourceOutputTransformer(EvalComponent, abc.ABC): ...
