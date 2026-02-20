from typing import Final, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "N_A",
    "Avogadro",
    "Boltzmann",
    "Btu",
    "Btu_IT",
    "Btu_th",
    "G",
    "Julian_year",
    "Planck",
    "R",
    "Rydberg",
    "Stefan_Boltzmann",
    "Wien",
    "acre",
    "alpha",
    "angstrom",
    "arcmin",
    "arcminute",
    "arcsec",
    "arcsecond",
    "astronomical_unit",
    "atm",
    "atmosphere",
    "atomic_mass",
    "atto",
    "au",
    "bar",
    "barrel",
    "bbl",
    "blob",
    "c",
    "calorie",
    "calorie_IT",
    "calorie_th",
    "carat",
    "centi",
    "convert_temperature",
    "day",
    "deci",
    "degree",
    "degree_Fahrenheit",
    "deka",
    "dyn",
    "dyne",
    "e",
    "eV",
    "electron_mass",
    "electron_volt",
    "elementary_charge",
    "epsilon_0",
    "erg",
    "exa",
    "exbi",
    "femto",
    "fermi",
    "fine_structure",
    "fluid_ounce",
    "fluid_ounce_US",
    "fluid_ounce_imp",
    "foot",
    "g",
    "gallon",
    "gallon_US",
    "gallon_imp",
    "gas_constant",
    "gibi",
    "giga",
    "golden",
    "golden_ratio",
    "grain",
    "gram",
    "gravitational_constant",
    "h",
    "hbar",
    "hectare",
    "hecto",
    "horsepower",
    "hour",
    "hp",
    "inch",
    "k",
    "kgf",
    "kibi",
    "kilo",
    "kilogram_force",
    "kmh",
    "knot",
    "lambda2nu",
    "lb",
    "lbf",
    "light_year",
    "liter",
    "litre",
    "long_ton",
    "m_e",
    "m_n",
    "m_p",
    "m_u",
    "mach",
    "mebi",
    "mega",
    "metric_ton",
    "micro",
    "micron",
    "mil",
    "mile",
    "milli",
    "minute",
    "mmHg",
    "mph",
    "mu_0",
    "nano",
    "nautical_mile",
    "neutron_mass",
    "nu2lambda",
    "ounce",
    "oz",
    "parsec",
    "pebi",
    "peta",
    "pi",
    "pico",
    "point",
    "pound",
    "pound_force",
    "proton_mass",
    "psi",
    "pt",
    "quecto",
    "quetta",
    "ronna",
    "ronto",
    "short_ton",
    "sigma",
    "slinch",
    "slug",
    "speed_of_light",
    "speed_of_sound",
    "stone",
    "survey_foot",
    "survey_mile",
    "tebi",
    "tera",
    "ton_TNT",
    "torr",
    "troy_ounce",
    "troy_pound",
    "u",
    "week",
    "yard",
    "year",
    "yobi",
    "yocto",
    "yotta",
    "zebi",
    "zepto",
    "zero_Celsius",
    "zetta",
]

_InexactArrayLikeT = TypeVar("_InexactArrayLikeT", bound=npc.inexact | onp.ArrayND[npc.inexact])

_TempScaleC: TypeAlias = Literal["Celsius", "celsius", "C", "c"]
_TempScaleK: TypeAlias = Literal["Kelvin", "kelvin", "K", "k"]
_TempScaleF: TypeAlias = Literal["Fahrenheit", "fahrenheit", "F", "f"]
_TempScaleR: TypeAlias = Literal["Rankine", "rankine", "R", "r"]
_TempScale: TypeAlias = Literal[_TempScaleC, _TempScaleK, _TempScaleF, _TempScaleR]

###

# mathematical constants
pi: Final = 3.141592653589793
golden: Final = 1.618033988749895
golden_ratio: Final = 1.618033988749895

# SI prefixes
quetta: Final = 1e30
ronna: Final = 1e27
yotta: Final = 1e24
zetta: Final = 1e21
exa: Final = 1e18
peta: Final = 1e15
tera: Final = 1e12
giga: Final = 1e9
mega: Final = 1e6
kilo: Final = 1e3
hecto: Final = 1e2
deka: Final = 1e1
deci: Final = 1e-1
centi: Final = 1e-2
milli: Final = 1e-3
micro: Final = 1e-6
nano: Final = 1e-9
pico: Final = 1e-12
femto: Final = 1e-15
atto: Final = 1e-18
zepto: Final = 1e-21
yocto: Final = 1e-24
ronto: Final = 1e-27
quecto: Final = 1e-30

# binary prefixes
kibi: Final[int] = ...  # 0x400
mebi: Final[int] = ...  # 0x10_0000
gibi: Final[int] = ...  # 0x4000_0000
tebi: Final[int] = ...  # 0x100_0000_0000
pebi: Final[int] = ...  # 0x4_0000_0000_0000
exbi: Final[int] = ...  # 0x1000_0000_0000_0000
zebi: Final[int] = ...  # 0x40_0000_0000_0000_0000
yobi: Final[int] = ...  # 0x1_0000_0000_0000_0000_0000

# physical constants
c: Final = 299792458.0
speed_of_light: Final = 299792458.0
mu_0: Final = 1.25663706212e-6
epsilon_0: Final = 8.8541878128e-12
h: Final = 6.62607015e-34
Planck: Final = 6.62607015e-34
hbar: Final = 1.0545718176461565e-34
G: Final = 6.6743e-11
gravitational_constant: Final = 6.6743e-11
g: Final = 9.80665
e: Final = 1.602176634e-19
elementary_charge: Final = 1.602176634e-19
R: Final = 8.314462618
gas_constant: Final = 8.314462618
alpha: Final = 0.0072973525693
fine_structure: Final = 0.0072973525693
N_A: Final = 6.02214076e23
Avogadro: Final = 6.02214076e23
k: Final = 1.380649e-23
Boltzmann: Final = 1.380649e-23
sigma: Final = 5.670374419e-8
Stefan_Boltzmann: Final = 5.670374419e-8
Wien: Final = 0.002897771955
Rydberg: Final = 10973731.56816

# mass in kg
gram: Final = 1e-3
metric_ton: Final = 1e3
# legacy mass units
grain: Final = 6.479891e-5
lb: Final = 0.45359237
pound: Final = 0.45359237
blob: Final = 175.12683524647636
slinch: Final = 175.12683524647636
slug: Final = 14.593902937206364
oz: Final = 0.028349523125
ounce: Final = 0.028349523124
stone: Final = 6.35029318
long_ton: Final = 1.0160469088e3
short_ton: Final = 907.18474
troy_ounce: Final = 0.0311034768
troy_pound: Final = 0.3732417216
carat: Final = 2e-4
# fundamental masses
m_e: Final = 9.1093837015e-31
electron_mass: Final = 9.1093837015e-31
m_p: Final = 1.67262192369e-27
proton_mass: Final = 1.67262192369e-27
m_n: Final = 1.67492749804e-27
neutron_mass: Final = 1.67492749804e-27
m_u: Final = 1.6605390666e-27
u: Final = 1.6605390666e-27
atomic_mass: Final = 1.6605390666e-27

# angle in rad
degree: Final = 0.017453292519943295
arcmin: Final = 2.908882086657216e-4
arcminute: Final = 2.908882086657216e-4
arcsec: Final = 4.84813681109536e-6
arcsecond: Final = 4.84813681109536e-6

# time in second
minute: Final = 60.0
hour: Final = 3_600.0
day: Final = 86_400.0
week: Final = 604_800.0
year: Final = 31_536_000.0
Julian_year: Final = 31_557_600.0

# length in meter
inch: Final = 0.0254
foot: Final = 0.3048
yard: Final = 0.9144
mile: Final = 1_609.344
mil: Final = 2.54e-5
pt: Final = 3.527777777777778e-4
point: Final = 3.527777777777778e-4
survey_foot: Final = 0.3048006096012192
survey_mile: Final = 1_609.3472186944373
nautical_mile: Final = 1_852.0
fermi: Final = 1e-15
angstrom: Final = 1e-10
micron: Final = 1e-6
au: Final = 149_597_870_700.0
astronomical_unit: Final = 149_597_870_700.0
light_year: Final = 9.460730472580800e15
parsec: Final = 3.085677581491367e16

# pressure in pascal
bar: Final = 100_000.0
atm: Final = 101_325.0
atmosphere: Final = 101_325.0
torr: Final = 133.32236842105263
mmHg: Final = 133.32236842105263
psi: Final = 6_894.757293168361

# area in meter**2
hectare: Final = 10_000.0
acre: Final = 4_046.8564224

# volume in meter**3
litre: Final = 1e-3
liter: Final = 1e-3
gallon: Final = 3.785411784e-3
gallon_US: Final = 3.785411784e-3
fluid_ounce: Final = 2.95735295625e-5
fluid_ounce_US: Final = 2.95735295625e-5
bbl: Final = 0.158987294928
barrel: Final = 0.158987294928
gallon_imp: Final = 4.54609e-3
fluid_ounce_imp: Final = 2.84130625e-5

# speed in meter per second
kmh: Final = 0.27777777777777777
mph: Final = 0.44704
mach: Final = 340.5
speed_of_sound: Final = 340.5
knot: Final = 0.5144444444444444

# temperature in kelvin
zero_Celsius: Final = 273.15
degree_Fahrenheit: Final = 0.55555555555555555

# energy in joule
eV: Final = 1.602176634e-19
electron_volt: Final = 1.602176634e-19
calorie: Final = 4.184
calorie_th: Final = 4.184
calorie_IT: Final = 4.1868
erg: Final = 1e-7
Btu_th: Final = 1_054.3502644888888
Btu: Final = 1_055.05585262
Btu_IT: Final = 1_055.05585262
ton_TNT: Final = 4.184e9

# power in watt
hp: Final = 745.6998715822701
horsepower: Final = 745.6998715822701

# force in newton
dyn: Final = 1e-5
dyne: Final = 1e-5
lbf: Final = 4.4482216152605
pound_force: Final = 4.4482216152605
kgf: Final = 9.80665
kilogram_force: Final = 9.80665

@overload
def convert_temperature(val: _InexactArrayLikeT, old_scale: _TempScale, new_scale: _TempScale) -> _InexactArrayLikeT: ...
@overload
def convert_temperature(val: onp.ToInt | onp.ToJustFloat64, old_scale: _TempScale, new_scale: _TempScale) -> np.float64: ...
@overload
def convert_temperature(
    val: onp.ToIntND | onp.ToJustFloat64_ND, old_scale: _TempScale, new_scale: _TempScale
) -> onp.ArrayND[np.float64]: ...
@overload
def convert_temperature(val: onp.ToFloatND, old_scale: _TempScale, new_scale: _TempScale) -> onp.ArrayND[npc.floating]: ...

#
@overload
def lambda2nu(lambda_: _InexactArrayLikeT) -> _InexactArrayLikeT: ...
@overload
def lambda2nu(lambda_: onp.ToInt | onp.ToJustFloat64) -> np.float64: ...
@overload
def lambda2nu(lambda_: onp.ToIntND | onp.ToJustFloat64_ND) -> onp.ArrayND[np.float64]: ...
@overload
def lambda2nu(lambda_: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...

#
@overload
def nu2lambda(nu: _InexactArrayLikeT) -> _InexactArrayLikeT: ...
@overload
def nu2lambda(nu: onp.ToInt | onp.ToJustFloat64) -> np.float64: ...
@overload
def nu2lambda(nu: onp.ToIntND | onp.ToJustFloat64_ND) -> onp.ArrayND[np.float64]: ...
@overload
def nu2lambda(nu: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...
