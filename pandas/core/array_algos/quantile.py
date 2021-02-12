import numpy as np

from pandas._libs import lib

from pandas.core.dtypes.common import is_list_like

from pandas.core.nanops import nanpercentile


def quantile_with_mask(
    values: np.ndarray,
    mask: np.ndarray,
    fill_value,
    qs,
    interpolation: str,
    axis: int,
) -> np.ndarray:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray
        For ExtensionArray, this is _values_for_factorize()[0]
    mask : np.ndarray[bool]
        mask = isna(values)
        For ExtensionArray, this is computed before calling _value_for_factorize
    fill_value : Scalar
        The value to interpret fill NA entries with
        For ExtensionArray, this is _values_for_factorize()[1]
    qs : a scalar or list of the quantiles to be computed
    interpolation : str
        Type of interpolation
    axis : int
        Axis along which to compute quantiles.

    Returns
    -------
    np.ndarray

    Notes
    -----
    Assumes values is already 2D.  For ExtensionArray this means np.atleast_2d
    has been called on _values_for_factorize()[0]
    """
    is_empty = values.shape[axis] == 0
    orig_scalar = not is_list_like(qs)
    if orig_scalar:
        # make list-like, unpack later
        qs = [qs]

    if is_empty:
        # create the array of na_values
        # 2d len(values) * len(qs)
        flat = np.array([fill_value] * len(qs))
        result = np.repeat(flat, len(values)).reshape(len(values), len(qs))
    else:
        # asarray needed for Sparse, see GH#24600
        result = nanpercentile(
            values,
            np.array(qs) * 100,
            axis=axis,
            na_value=fill_value,
            mask=mask,
            ndim=values.ndim,
            interpolation=interpolation,
        )

        result = np.array(result, copy=False)
        result = result.T

    if orig_scalar:
        assert result.shape[-1] == 1, result.shape
        result = result[..., 0]
        result = lib.item_from_zerodim(result)

    return result
