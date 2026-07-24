# ruff: noqa: PYI042
from typing import Literal as L  # noqa: N817

###

# boolean

type b1_name = L["bool", "bool_"]  # 'bool0' was removed in NumPy 2.0
type b1_char = L["b1", "|b1", "?"]
type b1_code = L[b1_name, b1_char]

# integer

type i1_name = L["int8", "byte"]
type i1_char = L["i1", "|i1", "b"]
type i1_code = L[i1_name, i1_char]

type u1_name = L["uint8", "ubyte"]
type u1_char = L["u1", "|u1", "B"]
type u1_code = L[u1_name, u1_char]

type i2_name = L["int16", "short"]
type i2_char = L["i2", "<i2", ">i2", "h"]
type i2_code = L[i2_name, i2_char]

type u2_name = L["uint16", "ushort"]
type u2_char = L["u2", "<u2", ">u2", "H"]
type u2_code = L[u2_name, u2_char]

type i4_name = L["int32", "intc"]
type i4_char = L["i4", "<i4", ">i4", "i"]
type i4_code = L[i4_name, i4_char]

type u4_name = L["uint32", "uintc"]
type u4_char = L["u4", "<u4", ">u4", "I"]
type u4_code = L[u4_name, u4_char]

type i8_name = L["int64", "longlong"]
type i8_char = L["i8", "<i8", ">i8", "q"]
type i8_code = L[i8_name, i8_char]

type u8_name = L["uint64", "ulonglong"]
type u8_char = L["u8", "<u8", ">u8", "Q"]
type u8_code = L[u8_name, u8_char]

type l_char = L["l", "<l", ">l"]
type L_char = L["L", "<L", ">L"]

type p_char = L["p", "<p", ">p"]
type P_char = L["P", "<P", ">P"]

# numpy>=2 only
type n_char = L["n", "<n", ">n"]
type N_char = L["N", "<N", ">N"]


type i0_name = L["int_", "int", "intp"]
type i0_char = n_char
type i0_code = L[i0_name, n_char]  # using `i0_char` here confuses ty

type u0_char = N_char
type u0_name = L["uint", "uintp"]
type u0_code = L[u0_name, N_char]  # using `u0_char` here confuses ty

type i__name = i0_name
type i__char = i0_char
type i__code = i0_code

type u__name = u0_name
type u__char = u0_char
type u__code = u0_code

type l_name = L["long"]
type l_code = L[l_name, l_char]

type L_name = L["ulong"]
type L_code = L[L_name, L_char]


# float

type f2_name = L["float16", "half"]
type f2_char = L["f2", "<f2", ">f2", "e"]
type f2_code = L[f2_name, f2_char]

type f4_name = L["float32", "single"]
type f4_char = L["f4", "<f4", ">f4", "f"]
type f4_code = L[f4_name, f4_char]

type f8_name = L["float64", "double", "float"]
type f8_char = L["f8", "<f8", ">f8", "d"]
type f8_code = L[f8_name, f8_char]

type f12_name = L["float96"]
type f12_char = L["f12", "<f12", ">f12"]
type f12_code = L[f12_name, f12_char]

type f16_name = L["float128"]
type f16_char = L["f16", "<f16", ">f16"]
type f16_code = L[f16_name, f16_char]

type g_name = L[f12_name, f16_name, "longdouble"]
type g_char = L[f12_char, f16_char, "g"]
type g_code = L[g_name, g_char]

# complex

type c8_name = L["complex64", "csingle"]
type c8_char = L["c8", "<c8", ">c8", "F"]
type c8_code = L[c8_name, c8_char]

type c16_name = L["complex", "complex128", "cdouble"]
type c16_char = L["c16", "<c16", ">c16", "D"]
type c16_code = L[c16_name, c16_char]

type c24_name = L["complex192"]
type c24_char = L["c24", "<c24", ">c24"]
type c24_code = L[c24_name, c24_char]

type c32_name = L["complex256"]
type c32_char = L["c32", "<c32", ">c32"]
type c32_code = L[c32_name, c32_char]

type G_name = L[c24_name, c32_name, "clongdouble"]
type G_char = L[c24_code, c32_code, "G"]
type G_code = L[G_name, G_char]

# object

type O_name = L["object_", "object"]
type O_char = L["O", "|O"]
type O_code = L[O_name, O_char]

# bytes_

# NOTE: this only includes 0-length bytes
type S0_name = L["bytes_", "bytes"]
type S0_char = L["S0", "|S0", "<S0", ">S0", "S"]
type S0_code = L[S0_name, S0_char]

type S1_name = L["bytes8"]
type S1_char = L["S1", "|S1", "<S1", ">S1", "c"]
type S1_code = L[S1_name, S1_char]

# str_

# NOTE: this only includes 0-length strings
type U0_name = L["str_", "str", "unicode"]
type U0_char = L["U0", "|U0", "<U0", ">U0", "U"]
type U0_code = L[U0_name, U0_char]

# void

# NOTE: this only includes "len-0 bytes void"
type V0_name = L["void"]  # 'void0' was removed in NumPy 2.0
type V0_char = L["V0", "|V0", "V"]
type V0_code = L[V0_name, V0_char]

# datetime64

type M_name = L[
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
type M_char = L[
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
type M_code = L[M_name, M_char]

# timedelta64

type m_name = L[
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
type m_char = L[
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
type m_code = L[m_name, m_char]

# stringv (or whatever we're gonna call the `StringDType().type` scalar type)

type T_name = L["StringDType128"]
type T_char = L["T"]
type T_code = T_char  # not yet

# abstract

type ix_name = L[
    "int8", "int16", "int32", "int64",
    "byte", "short", "intc", "long", "longlong",
    "int_", "intp",
]  # fmt: skip
type ix_char = L[i1_char, i2_char, i4_char, i8_char, l_char, p_char, n_char]
type ix_code = L[ix_name, ix_char]


type ux_name = L[
    "uint8", "uint16", "uint32", "uint64",
    "ubyte", "ushort", "uintc", "ulong", "ulonglong",
    "uint", "uintp",
]  # fmt: skip
type ux_char = L[u1_char, u2_char, u4_char, u8_char, L_char, P_char, N_char]
type ux_code = L[ux_name, ux_char]

type fx_code = L[f2_code, f4_code, f8_code, g_code]
type cx_code = L[c8_code, c16_code, G_code]

type iu_code = L[ux_code, ix_code]
type fc_code = L[fx_code, cx_code]
type iufc_code = L[iu_code, fc_code]

type SU_code = L[S0_code, S1_code, U0_code]
type SUV_code = L[SU_code, V0_code]
