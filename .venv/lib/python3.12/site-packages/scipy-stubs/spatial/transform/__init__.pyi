from . import rotation as rotation  # deprecated namespace, to be removed in v2.0.0
from ._rigid_transform import RigidTransform
from ._rotation import Rotation, Slerp
from ._rotation_spline import RotationSpline

__all__ = ["RigidTransform", "Rotation", "RotationSpline", "Slerp"]
