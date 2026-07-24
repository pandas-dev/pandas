import cmath


def div_usecase(x, y):
    return x / y


def real_usecase(x):
    return x.real

def imag_usecase(x):
    return x.imag

def conjugate_usecase(x):
    return x.conjugate()


def acos_usecase(x):
    return cmath.acos(x)

def cos_usecase(x):
    return cmath.cos(x)

def asin_usecase(x):
    return cmath.asin(x)

def sin_usecase(x):
    return cmath.sin(x)

def atan_usecase(x):
    return cmath.atan(x)

def tan_usecase(x):
    return cmath.tan(x)

def acosh_usecase(x):
    return cmath.acosh(x)

def cosh_usecase(x):
    return cmath.cosh(x)

def asinh_usecase(x):
    return cmath.asinh(x)

def sinh_usecase(x):
    return cmath.sinh(x)

def atanh_usecase(x):
    return cmath.atanh(x)

def tanh_usecase(x):
    return cmath.tanh(x)

def exp_usecase(x):
    return cmath.exp(x)

def isfinite_usecase(x):
    return cmath.isfinite(x)

def isinf_usecase(x):
    return cmath.isinf(x)

def isnan_usecase(x):
    return cmath.isnan(x)

def log_usecase(x):
    return cmath.log(x)

def log_base_usecase(x, base):
    return cmath.log(x, base)

def log10_usecase(x):
    return cmath.log10(x)

def phase_usecase(x):
    return cmath.phase(x)

def polar_usecase(x):
    return cmath.polar(x)

_two = 2.0

def polar_as_complex_usecase(x):
    # HACK: clear errno by invoking float.__pow__
    # (workaround for http://bugs.python.org/issue24489)
    _two ** _two
    return complex(*cmath.polar(x))

def rect_usecase(r, phi):
    return cmath.rect(r, phi)

def sqrt_usecase(x):
    return cmath.sqrt(x)
