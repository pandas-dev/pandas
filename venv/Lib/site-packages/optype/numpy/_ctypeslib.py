import ctypes as ct

__all__ = "CScalar", "CType"


_, CType, CScalar, *_ = ct.c_int.mro()
