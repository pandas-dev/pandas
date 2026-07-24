from typing import ClassVar

from .geometry import Point2D, PositiveSize2D, Transform2D

class XDRPoint2D(Point2D):
    namespace: ClassVar[None]  # type: ignore[assignment]
    # Same as parent
    # x = Point2D.x
    # y = Point2D.y

class XDRPositiveSize2D(PositiveSize2D):
    namespace: ClassVar[None]  # type: ignore[assignment]
    # Same as parent
    # cx = PositiveSize2D.cx
    # cy = PositiveSize2D.cy

class XDRTransform2D(Transform2D):
    namespace: ClassVar[None]  # type: ignore[assignment]
    # Same as parent
    # rot = Transform2D.rot
    # flipH = Transform2D.flipH
    # flipV = Transform2D.flipV
    # off = Transform2D.off
    # ext = Transform2D.ext
    # chOff = Transform2D.chOff
    # chExt = Transform2D.chExt
