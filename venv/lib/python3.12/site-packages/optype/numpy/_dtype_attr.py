# ruff: noqa: PYI042
from typing import Literal as L, TypeAlias as Alias  # noqa: N817

from numpy_typing_compat import NUMPY_GE_2_0

###

# boolean

b1_name: Alias = L["bool", "bool_"]  # 'bool0' was removed in NumPy 2.0
b1_char: Alias = L["b1", "|b1", "?"]
b1_code: Alias = L[b1_name, b1_char]

# integer

i1_name: Alias = L["int8", "byte"]
i1_char: Alias = L["i1", "|i1", "b"]
i1_code: Alias = L[i1_name, i1_char]

u1_name: Alias = L["uint8", "ubyte"]
u1_char: Alias = L["u1", "|u1", "B"]
u1_code: Alias = L[u1_name, u1_char]

i2_name: Alias = L["int16", "short"]
i2_char: Alias = L["i2", "<i2", ">i2", "h"]
i2_code: Alias = L[i2_name, i2_char]

u2_name: Alias = L["uint16", "ushort"]
u2_char: Alias = L["u2", "<u2", ">u2", "H"]
u2_code: Alias = L[u2_name, u2_char]

i4_name: Alias = L["int32", "intc"]
i4_char: Alias = L["i4", "<i4", ">i4", "i"]
i4_code: Alias = L[i4_name, i4_char]

u4_name: Alias = L["uint32", "uintc"]
u4_char: Alias = L["u4", "<u4", ">u4", "I"]
u4_code: Alias = L[u4_name, u4_char]

i8_name: Alias = L["int64", "longlong"]
i8_char: Alias = L["i8", "<i8", ">i8", "q"]
i8_code: Alias = L[i8_name, i8_char]

u8_name: Alias = L["uint64", "ulonglong"]
u8_char: Alias = L["u8", "<u8", ">u8", "Q"]
u8_code: Alias = L[u8_name, u8_char]

l_char: Alias = L["l", "<l", ">l"]
L_char: Alias = L["L", "<L", ">L"]

p_char: Alias = L["p", "<p", ">p"]
P_char: Alias = L["P", "<P", ">P"]

# numpy>=2 only
n_char: Alias = L["n", "<n", ">n"]
N_char: Alias = L["N", "<N", ">N"]


# mypy: disable-error-code="no-redef"
# pyright: reportRedeclaration=false

if NUMPY_GE_2_0:
    i0_name: Alias = L["int_", "int", "intp"]
    i0_char: Alias = n_char
    i0_code: Alias = L[i0_name, i0_char]

    u0_char: Alias = N_char
    u0_name: Alias = L["uint", "uintp"]
    u0_code: Alias = L[u0_name, u0_char]

    i__name: Alias = i0_name
    i__char: Alias = i0_char
    i__code: Alias = i0_code

    u__name: Alias = u0_name
    u__char: Alias = u0_char
    u__code: Alias = u0_code

    l_name: Alias = L["long"]
    l_code: Alias = L[l_name, l_char]

    L_name: Alias = L["ulong"]
    L_code: Alias = L[L_name, L_char]
else:
    i0_name: Alias = L["intp"]  # type: ignore[misc]
    i0_char: Alias = p_char  # type: ignore[misc]
    i0_code: Alias = L[i0_name, i0_char]  # type: ignore[misc]

    u0_name: Alias = L["uintp"]  # type: ignore[misc]
    u0_char: Alias = P_char  # type: ignore[misc]
    u0_code: Alias = L[u0_name, u0_char]  # type: ignore[misc]

    i__name: Alias = L["int_", "int", "long"]  # type: ignore[misc]
    i__char: Alias = L[l_char]  # type: ignore[misc]
    i__code: Alias = L[i__name, i__char]  # type: ignore[misc]

    u__name: Alias = L["uint", "ulong"]  # type: ignore[misc]
    u__char: Alias = L[L_char]  # type: ignore[misc]
    u__code: Alias = L[u__name, u__char]  # type: ignore[misc]

    l_name: Alias = i__name  # type: ignore[misc]
    l_code: Alias = i__code  # type: ignore[misc]

    L_name: Alias = u__name  # type: ignore[misc]
    L_code: Alias = u__code  # type: ignore[misc]


# float

f2_name: Alias = L["float16", "half"]
f2_char: Alias = L["f2", "<f2", ">f2", "e"]
f2_code: Alias = L[f2_name, f2_char]

f4_name: Alias = L["float32", "single"]
f4_char: Alias = L["f4", "<f4", ">f4", "f"]
f4_code: Alias = L[f4_name, f4_char]

f8_name: Alias = L["float64", "double", "float"]
f8_char: Alias = L["f8", "<f8", ">f8", "d"]
f8_code: Alias = L[f8_name, f8_char]

f12_name: Alias = L["float96"]
f12_char: Alias = L["f12", "<f12", ">f12"]
f12_code: Alias = L[f12_name, f12_char]

f16_name: Alias = L["float128"]
f16_char: Alias = L["f16", "<f16", ">f16"]
f16_code: Alias = L[f16_name, f16_char]

g_name: Alias = L[f12_name, f16_name, "longdouble"]
g_char: Alias = L[f12_char, f16_char, "g"]
g_code: Alias = L[g_name, g_char]

# complex

c8_name: Alias = L["complex64", "csingle"]
c8_char: Alias = L["c8", "<c8", ">c8", "F"]
c8_code: Alias = L[c8_name, c8_char]

c16_name: Alias = L["complex", "complex128", "cdouble"]
c16_char: Alias = L["c16", "<c16", ">c16", "D"]
c16_code: Alias = L[c16_name, c16_char]

c24_name: Alias = L["complex192"]
c24_char: Alias = L["c24", "<c24", ">c24"]
c24_code: Alias = L[c24_name, c24_char]

c32_name: Alias = L["complex256"]
c32_char: Alias = L["c32", "<c32", ">c32"]
c32_code: Alias = L[c32_name, c32_char]

G_name: Alias = L[c24_name, c32_name, "clongdouble"]
G_char: Alias = L[c24_code, c32_code, "G"]
G_code: Alias = L[G_name, G_char]

# object

O_name: Alias = L["object_", "object"]
O_char: Alias = L["O", "|O"]
O_code: Alias = L[O_name, O_char]

# bytes_

# NOTE: this only includes 0-length bytes
S0_name: Alias = L["bytes_", "bytes"]
S0_char: Alias = L["S0", "|S0", "<S0", ">S0", "S"]
S0_code: Alias = L[S0_name, S0_char]

S1_name: Alias = L["bytes8"]
S1_char: Alias = L["S1", "|S1", "<S1", ">S1", "c"]
S1_code: Alias = L[S1_name, S1_char]

# str_

# NOTE: this only includes 0-length strings
U0_name: Alias = L["str_", "str", "unicode"]
U0_char: Alias = L["U0", "|U0", "<U0", ">U0", "U"]
U0_code: Alias = L[U0_name, U0_char]

# void

# NOTE: this only includes "len-0 bytes void"
V0_name: Alias = L["void"]  # 'void0' was removed in NumPy 2.0
V0_char: Alias = L["V0", "|V0", "V"]
V0_code: Alias = L[V0_name, V0_char]

# datetime64

M_name: Alias = L[
    "datetime64",
    "datetime64[as]",
    "datetime64[fs]",
    "datetime64[ps]",
    "datetime64[ns]",
    "datetime64[us]",
    "datetime64[ms]",
    "datetime64[s]",
    "datetime64[m]",
    "datetime64[h]",
    "datetime64[D]",
    "datetime64[W]",
    "datetime64[M]",
    "datetime64[Y]",
]
M_char: Alias = L[
    "M8", "<M8", ">M8", "M",
    "M8[as]", "<M8[as]", ">M8[as]",
    "M8[fs]", "<M8[fs]", ">M8[fs]",
    "M8[ps]", "<M8[ps]", ">M8[ps]",
    "M8[ns]", "<M8[ns]", ">M8[ns]",
    "M8[us]", "<M8[us]", ">M8[us]",
    "M8[s]", "<M8[s]", ">M8[s]",
    "M8[m]", "<M8[m]", ">M8[m]",
    "M8[h]", "<M8[h]", ">M8[h]",
    "M8[D]", "<M8[D]", ">M8[D]",
    "M8[W]", "<M8[W]", ">M8[W]",
    "M8[M]", "<M8[M]", ">M8[M]",
    "M8[Y]", "<M8[Y]", ">M8[Y]",
]  # fmt: skip
M_code: Alias = L[M_name, M_char]

# timedelta64

m_name: Alias = L[
    "timedelta64",
    "timedelta64[as]",
    "timedelta64[fs]",
    "timedelta64[ps]",
    "timedelta64[ns]",
    "timedelta64[us]",
    "timedelta64[ms]",
    "timedelta64[s]",
    "timedelta64[m]",
    "timedelta64[h]",
    "timedelta64[D]",
    "timedelta64[W]",
    "timedelta64[M]",
    "timedelta64[Y]",
]
m_char: Alias = L[
    "m8", "<m8", ">m8", "m",
    "m8[as]", "<m8[as]", ">m8[as]",
    "m8[fs]", "<m8[fs]", ">m8[fs]",
    "m8[ps]", "<m8[ps]", ">m8[ps]",
    "m8[ns]", "<m8[ns]", ">m8[ns]",
    "m8[us]", "<m8[us]", ">m8[us]",
    "m8[s]", "<m8[s]", ">m8[s]",
    "m8[m]", "<m8[m]", ">m8[m]",
    "m8[h]", "<m8[h]", ">m8[h]",
    "m8[D]", "<m8[D]", ">m8[D]",
    "m8[W]", "<m8[W]", ">m8[W]",
    "m8[M]", "<m8[M]", ">m8[M]",
    "m8[Y]", "<m8[Y]", ">m8[Y]",
]  # fmt: skip
m_code: Alias = L[m_name, m_char]

# stringv (or whatever we're gonna call the `StringDType().type` scalar type)

T_name: Alias = L["StringDType128"]
T_char: Alias = L["T"]
T_code: Alias = T_char  # not yet

# abstract

ix_name: Alias = L[
    "int8", "int16", "int32", "int64",
    "byte", "short", "intc", "long", "longlong",
    "int_", "intp",
]  # fmt: skip
ix_char: Alias = L[i1_char, i2_char, i4_char, i8_char, l_char, p_char, n_char]
ix_code: Alias = L[ix_name, ix_char]


ux_name: Alias = L[
    "uint8", "uint16", "uint32", "uint64",
    "ubyte", "ushort", "uintc", "ulong", "ulonglong",
    "uint", "uintp",
]  # fmt: skip
ux_char: Alias = L[u1_char, u2_char, u4_char, u8_char, L_char, P_char, N_char]
ux_code: Alias = L[ux_name, ux_char]

fx_code: Alias = L[f2_code, f4_code, f8_code, g_code]
cx_code: Alias = L[c8_code, c16_code, G_code]

iu_code: Alias = L[ux_code, ix_code]
fc_code: Alias = L[fx_code, cx_code]
iufc_code: Alias = L[iu_code, fc_code]

SU_code: Alias = L[S0_code, S1_code, U0_code]
SUV_code: Alias = L[SU_code, V0_code]
