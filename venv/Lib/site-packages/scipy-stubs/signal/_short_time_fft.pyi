import types
from collections.abc import Callable
from typing import Any, Final, Generic, Literal, Never, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from .windows._windows import _ToWindow

__all__ = ["ShortTimeFFT", "closest_STFT_dual_window"]

###

_InexactT_co = TypeVar("_InexactT_co", bound=npc.inexact, default=Any, covariant=True)
# `float64` if `fft_mode in {"onesided", "onesided2X"}` otherwise `complex128`; only affects `istft`
_Inexact128T_co = TypeVar("_Inexact128T_co", bound=np.float64 | np.complex128, default=Any, covariant=True)

type _PadType = Literal["zeros", "edge", "even", "odd"]

type _FFTMode1 = Literal["onesided", "onesided2X"]
type _FFTMode2 = Literal["twosided", "centered"]
type _FFTMode = Literal[_FFTMode1, _FFTMode2]

type _ScaleTo = Literal["magnitude", "psd"]
type _Scaling = Literal[_ScaleTo, "unitary"]

type _Detr = (
    Literal["linear", "constant"]
    | Callable[[onp.ArrayND[np.float64]], onp.ToComplexND]
    | Callable[[onp.ArrayND[np.complex128]], onp.ToComplexND]
)

###

class ShortTimeFFT(Generic[_InexactT_co, _Inexact128T_co]):
    _win: onp.Array1D[_InexactT_co]
    _dual_win: onp.Array1D[_InexactT_co] | None = None
    _hop: Final[int]

    _fs: float
    _fft_mode: _FFTMode = "onesided"
    _mfft: int
    _scaling: _Scaling | None = None
    _phase_shift: int | None

    _fac_mag: float | None = None
    _fac_psd: float | None = None
    _lower_border_end: tuple[int, int] | None = None

    _cache_post_padding: tuple[int, tuple[int, int]] = ...
    _cache_upper_border_begin: tuple[int, tuple[int, int]] = ...
    _cache_t: tuple[tuple[int, int | None, int | None, int, float], onp.Array1D[np.float64]] = ...
    _cache_f: tuple[tuple[_FFTMode, int, float], onp.Array1D[np.float64]] = ...

    @classmethod
    def __class_getitem__(cls, arg: object, /) -> types.GenericAlias: ...

    #
    @overload
    @classmethod
    def from_dual[InexactT: npc.inexact](
        cls,
        dual_win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode1 = "onesided",
        mfft: int | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[InexactT, np.float64]: ...
    @overload
    @classmethod
    def from_dual[InexactT: npc.inexact](
        cls,
        dual_win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode2,
        mfft: int | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[InexactT, np.complex128]: ...

    #
    @overload
    @classmethod
    def from_window(
        cls,
        win_param: _ToWindow,
        fs: float,
        nperseg: int,
        noverlap: int,
        *,
        symmetric_win: bool = False,
        fft_mode: _FFTMode1 = "onesided",
        mfft: int | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[np.float64, np.float64]: ...
    @overload
    @classmethod
    def from_window(
        cls,
        win_param: _ToWindow,
        fs: float,
        nperseg: int,
        noverlap: int,
        *,
        symmetric_win: bool = False,
        fft_mode: _FFTMode2,
        mfft: int | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[np.float64, np.complex128]: ...

    #
    @overload
    @classmethod
    def from_win_equals_dual[InexactT: npc.inexact](
        cls,
        desired_win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode1 = "onesided",
        mfft: int | None = None,
        scale_to: _Scaling | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[InexactT, np.float64]: ...
    @overload
    @classmethod
    def from_win_equals_dual[InexactT: npc.inexact](
        cls,
        desired_win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode2,
        mfft: int | None = None,
        scale_to: _Scaling | None = None,
        phase_shift: int | None = 0,
    ) -> ShortTimeFFT[InexactT, np.complex128]: ...

    #
    @property
    def win(self, /) -> onp.Array1D[_InexactT_co]: ...
    @property
    def dual_win(self, /) -> onp.Array1D[_InexactT_co]: ...
    @property
    def hop(self, /) -> int: ...
    @property
    def invertible(self, /) -> bool: ...
    @property
    def fac_magnitude(self, /) -> float: ...
    @property
    def fac_psd(self, /) -> float: ...
    @property
    def m_num(self, /) -> int: ...
    @property
    def m_num_mid(self, /) -> int: ...
    @property
    def k_min(self, /) -> int: ...
    @property
    def p_min(self, /) -> int: ...
    @property
    def lower_border_end(self, /) -> tuple[int, int]: ...
    @property
    def delta_t(self, /) -> float: ...
    @property
    def delta_f(self, /) -> float: ...
    @property
    def f_pts(self, /) -> int: ...
    @property
    def f(self, /) -> onp.Array1D[np.float64]: ...
    @property
    def onesided_fft(self, /) -> bool: ...
    @property
    def scaling(self, /) -> _Scaling | None: ...

    #
    @property
    def T(self, /) -> float: ...
    @T.setter
    def T(self, /, v: float) -> None: ...

    #
    @property
    def fs(self, /) -> float: ...
    @fs.setter
    def fs(self, /, v: float) -> None: ...

    #
    @property
    def fft_mode(self, /) -> _FFTMode: ...
    @fft_mode.setter  # `fft_mode` affects the 2nd generic type (`_Inexact128T_co`), so changing it would be type-unsafe
    def fft_mode(self, /, t: Never) -> None: ...

    #
    @property
    def mfft(self, /) -> int: ...
    @mfft.setter
    def mfft(self, /, n_: int) -> None: ...

    #
    @property
    def phase_shift(self, /) -> int | None: ...
    @phase_shift.setter
    def phase_shift(self, /, v: int | None) -> None: ...

    #
    @property
    def _pre_padding(self, /) -> tuple[int, int]: ...
    def _post_padding(self, /, n: int) -> tuple[int, int]: ...

    #
    @overload
    def __init__[InexactT: npc.inexact](
        self: ShortTimeFFT[InexactT, np.float64],
        /,
        win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode1 = "onesided",
        mfft: int | None = None,
        dual_win: onp.ArrayND[InexactT] | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> None: ...
    @overload
    def __init__[InexactT: npc.inexact](
        self: ShortTimeFFT[InexactT, np.complex128],
        /,
        win: onp.ArrayND[InexactT],
        hop: int,
        fs: float,
        *,
        fft_mode: _FFTMode2,
        mfft: int | None = None,
        dual_win: onp.ArrayND[InexactT] | None = None,
        scale_to: _ScaleTo | None = None,
        phase_shift: int | None = 0,
    ) -> None: ...

    #
    def k_max(self, /, n: int) -> int: ...
    def p_max(self, /, n: int) -> int: ...
    def p_num(self, /, n: int) -> int: ...
    def nearest_k_p(self, /, k: int, left: bool = True) -> int: ...
    def upper_border_begin(self, /, n: int) -> tuple[int, int]: ...
    def p_range(self, /, n: int, p0: int | None = None, p1: int | None = None) -> tuple[int, int]: ...
    def t(self, /, n: int, p0: int | None = None, p1: int | None = None, k_offset: int = 0) -> onp.Array1D[np.float64]: ...
    def scale_to(self, /, scaling: _ScaleTo) -> None: ...

    #
    @overload
    def stft(
        self,
        /,
        x: onp.Array1D[npc.inexact],
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array2D[np.complex128]: ...
    @overload
    def stft(
        self,
        /,
        x: onp.Array2D[npc.inexact],
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array3D[np.complex128]: ...
    @overload
    def stft(
        self,
        /,
        x: onp.ArrayND[npc.inexact],
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.ArrayND[np.complex128]: ...

    #
    @overload
    def stft_detrend(
        self,
        /,
        x: onp.Array1D[npc.inexact],
        detr: _Detr | None,
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array2D[np.complex128]: ...
    @overload
    def stft_detrend(
        self,
        /,
        x: onp.Array2D[npc.inexact],
        detr: _Detr | None,
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array3D[np.complex128]: ...
    @overload
    def stft_detrend(
        self,
        /,
        x: onp.ArrayND[npc.inexact],
        detr: _Detr | None,
        p0: int | None = None,
        p1: int | None = None,
        *,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.ArrayND[np.complex128]: ...

    #
    @overload
    def spectrogram(
        self,
        /,
        x: onp.Array1D[npc.inexact],
        y: None = None,
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array2D[np.float64]: ...
    @overload
    def spectrogram(
        self,
        /,
        x: onp.Array1D[npc.inexact],
        y: onp.Array1D[npc.inexact],
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array2D[np.complex128]: ...
    @overload
    def spectrogram(
        self,
        /,
        x: onp.Array2D[npc.inexact],
        y: None = None,
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array3D[np.float64]: ...
    @overload
    def spectrogram(
        self,
        /,
        x: onp.Array2D[npc.inexact],
        y: onp.Array2D[npc.inexact],
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.Array3D[np.complex128]: ...
    @overload
    def spectrogram(
        self,
        /,
        x: onp.ArrayND[npc.inexact],
        y: None = None,
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def spectrogram(
        self,
        /,
        x: onp.ArrayND[npc.inexact],
        y: onp.ArrayND[npc.inexact],
        detr: _Detr | None = None,
        *,
        p0: int | None = None,
        p1: int | None = None,
        k_offset: int = 0,
        padding: _PadType = "zeros",
        axis: int = -1,
    ) -> onp.ArrayND[np.complex128]: ...

    #
    def istft(
        self, /, S: onp.ArrayND[npc.inexact], k0: int = 0, k1: int | None = None, *, f_axis: int = -2, t_axis: int = -1
    ) -> onp.ArrayND[_Inexact128T_co]: ...

    #
    def extent(
        self, /, n: int, axes_seq: Literal["tf", "ft"] = "tf", center_bins: bool = False
    ) -> tuple[float, float, float, float]: ...

#
def _calc_dual_canonical_window[InexactT: npc.inexact](win: onp.ArrayND[InexactT], hop: int) -> onp.Array1D[InexactT]: ...

#
@overload
def closest_STFT_dual_window[InexactT: npc.inexact](
    win: onp.ArrayND[InexactT], hop: int, desired_dual: onp.ArrayND[InexactT] | None = None, *, scaled: onp.ToTrue = True
) -> tuple[onp.Array1D[InexactT], InexactT]: ...
@overload
def closest_STFT_dual_window[InexactT: npc.inexact](
    win: onp.ArrayND[InexactT], hop: int, desired_dual: onp.ArrayND[InexactT] | None = None, *, scaled: onp.ToFalse
) -> tuple[onp.Array1D[InexactT], float]: ...
