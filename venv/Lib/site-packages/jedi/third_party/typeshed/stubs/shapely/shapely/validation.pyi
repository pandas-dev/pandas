from .geometry.base import BaseGeometry
from .lib import Geometry

__all__ = ["explain_validity", "make_valid"]

def explain_validity(ob: Geometry) -> str: ...
def make_valid(ob: Geometry) -> BaseGeometry: ...
