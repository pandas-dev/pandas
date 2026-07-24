from typing import Final, Literal, overload

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

###

type _TempScaleC = Literal["Celsius", "celsius", "C", "c"]
type _TempScaleK = Literal["Kelvin", "kelvin", "K", "k"]
type _TempScaleF = Literal["Fahrenheit", "fahrenheit", "F", "f"]
type _TempScaleR = Literal["Rankine", "rankine", "R", "r"]
type _TempScale = Literal[_TempScaleC, _TempScaleK, _TempScaleF, _TempScaleR]

###

# mathematical constants
pi: Final = 3.141592653589793
golden: Final = 1.618033988749895
golden_ratio: Final = golden

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
c: Final = 299_792_458.0
speed_of_light: Final = c
mu_0: Final = 1.25663706127e-6
epsilon_0: Final = 8.8541878188e-12
h: Final = 6.62607015e-34
Planck: Final = h
hbar: Final = 1.0545718176461565e-34
G: Final = 6.6743e-11
gravitational_constant: Final = G
g: Final = 9.80665
e: Final = 1.602176634e-19
elementary_charge: Final = e
R: Final = 8.31446261815324
gas_constant: Final = R
alpha: Final = 0.0072973525643
fine_structure: Final = alpha
N_A: Final = 6.02214076e23
Avogadro: Final = N_A
k: Final = 1.380649e-23
Boltzmann: Final = k
sigma: Final = 5.6703744191844314e-8
Stefan_Boltzmann: Final = sigma
Wien: Final = 0.0028977719551851727
Rydberg: Final = 10_973_731.568157

# mass in kg
gram: Final = 1e-3
metric_ton: Final = 1e3
# legacy mass units
grain: Final = 6.479891e-5
lb: Final = 0.45359236999999997
pound: Final = lb
blob: Final = 175.12683524647636
slinch: Final = blob
slug: Final = 14.593902937206364
oz: Final = 0.028349523124999998
ounce: Final = oz
stone: Final = 6.3502931799999995
long_ton: Final = 1.0160469088e3
short_ton: Final = 907.1847399999999
troy_ounce: Final = 0.031103476799999998
troy_pound: Final = 0.37324172159999996
carat: Final = 2e-4
# fundamental masses
m_e: Final = 9.1093837139e-31
electron_mass: Final = m_e
m_p: Final = 1.67262192595e-27
proton_mass: Final = m_p
m_n: Final = 1.67492750056e-27
neutron_mass: Final = m_n
m_u: Final = 1.66053906892e-27
u: Final = m_u
atomic_mass: Final = m_u

# angle in rad
degree: Final = 0.017453292519943295
arcmin: Final = 2.908882086657216e-4
arcminute: Final = arcmin
arcsec: Final = 4.84813681109536e-6
arcsecond: Final = arcsec

# time in second
minute: Final = 60.0
hour: Final = 3_600.0
day: Final = 86_400.0
week: Final = 604_800.0
year: Final = 31_536_000.0
Julian_year: Final = 31_557_600.0

# length in meter
inch: Final = 0.0254
foot: Final = 0.30479999999999996
yard: Final = 0.9143999999999999
mile: Final = 1_609.3439999999998
mil: Final = 2.5399999999999997e-5
pt: Final = 0.00035277777777777776
point: Final = pt
survey_foot: Final = 0.3048006096012192
survey_mile: Final = 1_609.3472186944373
nautical_mile: Final = 1_852.0
fermi: Final = 1e-15
angstrom: Final = 1e-10
micron: Final = 1e-6
au: Final = 149_597_870_700.0
astronomical_unit: Final = au
light_year: Final = 9.460730472580800e15
parsec: Final = 3.085677581491367e16

# pressure in pascal
bar: Final = 100_000.0
atm: Final = 101_325.0
atmosphere: Final = atm
torr: Final = 133.32236842105263
mmHg: Final = torr
psi: Final = 6_894.757293168361

# area in meter**2
hectare: Final = 10_000.0
acre: Final = 4_046.8564223999992

# volume in meter**3
litre: Final = 1e-3
liter: Final = litre
gallon: Final = 0.0037854117839999997
gallon_US: Final = gallon
fluid_ounce: Final = 2.9573529562499998e-5
fluid_ounce_US: Final = fluid_ounce
bbl: Final = 0.15898729492799998
barrel: Final = bbl
gallon_imp: Final = 4.54609e-3
fluid_ounce_imp: Final = 2.84130625e-5

# speed in meter per second
kmh: Final = 0.27777777777777777
mph: Final = 0.44703999999999994
mach: Final = 340.5
speed_of_sound: Final = mach
knot: Final = 0.5144444444444445

# temperature in kelvin
zero_Celsius: Final = 273.15
degree_Fahrenheit: Final = 0.55555555555555555

# energy in joule
eV: Final = 1.602176634e-19
electron_volt: Final = eV
calorie: Final = 4.184
calorie_th: Final = calorie
calorie_IT: Final = 4.1868
erg: Final = 1e-7
Btu_th: Final = 1_054.3502644888888
Btu: Final = 1_055.05585262
Btu_IT: Final = Btu
ton_TNT: Final = 4.184e9

# power in watt
hp: Final = 745.6998715822701
horsepower: Final = hp

# force in newton
dyn: Final = 1e-5
dyne: Final = dyn
lbf: Final = 4.4482216152605
pound_force: Final = lbf
kgf: Final = 9.80665
kilogram_force: Final = kgf

@overload
def convert_temperature[ArrayLikeT: npc.inexact | onp.ArrayND[npc.inexact]](
    val: ArrayLikeT, old_scale: _TempScale, new_scale: _TempScale
) -> ArrayLikeT: ...
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
def lambda2nu[ArrayLikeT: npc.inexact | onp.ArrayND[npc.inexact]](lambda_: ArrayLikeT) -> ArrayLikeT: ...
@overload
def lambda2nu(lambda_: onp.ToInt | onp.ToJustFloat64) -> np.float64: ...
@overload
def lambda2nu(lambda_: onp.ToIntND | onp.ToJustFloat64_ND) -> onp.ArrayND[np.float64]: ...
@overload
def lambda2nu(lambda_: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...

#
@overload
def nu2lambda[ArrayLikeT: npc.inexact | onp.ArrayND[npc.inexact]](nu: ArrayLikeT) -> ArrayLikeT: ...
@overload
def nu2lambda(nu: onp.ToInt | onp.ToJustFloat64) -> np.float64: ...
@overload
def nu2lambda(nu: onp.ToIntND | onp.ToJustFloat64_ND) -> onp.ArrayND[np.float64]: ...
@overload
def nu2lambda(nu: onp.ToFloatND) -> onp.ArrayND[npc.floating]: ...
