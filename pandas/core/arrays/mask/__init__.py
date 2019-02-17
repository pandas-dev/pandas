_MaskArrayType = None


def get_mask_array_type():
    """Set the mask array type to use, we need to do
    this after all modules are imported as the implementations
    e.g. pyarrow depend on pandas being importable
    """
    global _MaskArrayType

    if _MaskArrayType is not None:
        return _MaskArrayType

    # if ArrowBoolArray is available use it
    # otherwise use the NumpyMask
    try:
        from pandas.core.arrays.mask._pyarrow import ArrowMaskArray

        MaskArray = ArrowMaskArray

    except ImportError:
        from pandas.core.arrays.mask._numpy import NumpyMaskArray

        MaskArray = NumpyMaskArray

    _MaskArrayType = MaskArray
    return _MaskArrayType


__all__ = ['get_mask_array_type']
