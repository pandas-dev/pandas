#include <numpy/npy_cpu.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "numexpr_config.hpp" // isnan definitions

// Generic sign function
inline int signi(int x) {return (0 < x) - (x < 0);}
inline long signl(long x) {return (0 < x) - (x < 0);}
inline double sign(double x){
        // Floats: -1.0, 0.0, +1.0, NaN stays NaN
        if (isnand(x)) {return NAN;}
        if (x > 0) {return 1;}
        if (x < 0) {return -1;}
        return 0; // handles +0.0 and -0.0
    }
inline float signf(float x){
        // Floats: -1.0, 0.0, +1.0, NaN stays NaN
        if (isnanf_(x)) {return NAN;}
        if (x > 0) {return 1;}
        if (x < 0) {return -1;}
        return 0; // handles +0.0 and -0.0
    }

// round function for ints
inline int rinti(int x) {return x;}
inline long rintl(long x) {return x;}
// abs function for ints
inline int fabsi(int x) {return x<0 ? -x: x;}
inline long fabsl(long x) {return x<0 ? -x: x;}
// fmod function for ints
// TODO: Have to add FUNC_III, FUNC_LLL signatures to functions.hpp to enable these
// inline int fmodi(int x, int y) {return (int)fmodf((float)x, (float)y);}
// inline long fmodl(long x, long y)  {return (long)fmodf((long)x, (long)y);}

#ifdef USE_VML
// To match Numpy behaviour for NaNs
static void vsFmax_(MKL_INT n, const float* x1, const float* x2, float* dest)
{
    vsFmax(n, x1, x2, dest);
    MKL_INT j;
    for (j=0; j<n; j++) {
        if (isnanf_(x1[j]) | isnanf_(x2[j])){
            dest[j] = NAN;
        }
    };
};

static void vsFmin_(MKL_INT n, const float* x1, const float* x2, float* dest)
{
    vsFmin(n, x1, x2, dest);
    MKL_INT j;
    for (j=0; j<n; j++) {
        if (isnanf_(x1[j]) | isnanf_(x2[j])){
            dest[j] = NAN;
        }
    };
};
// To match Numpy behaviour for NaNs
static void vdFmax_(MKL_INT n, const double* x1, const double* x2, double* dest)
{
    vdFmax(n, x1, x2, dest);
    MKL_INT j;
    for (j=0; j<n; j++) {
        if (isnand(x1[j]) | isnand(x2[j])){
            dest[j] = NAN;
        }
    };
};

static void vdFmin_(MKL_INT n, const double* x1, const double* x2, double* dest)
{
    vdFmin(n, x1, x2, dest);
    MKL_INT j;
    for (j=0; j<n; j++) {
        if (isnand(x1[j]) | isnand(x2[j])){
            dest[j] = NAN;
        }
    };
};

static void viRint(MKL_INT n, const int* x, int* dest)
{
    memcpy(dest, x, n * sizeof(int)); // just copy x1 which is already int
};

static void vlRint(MKL_INT n, const long* x, long* dest)
{
    memcpy(dest, x, n * sizeof(long)); // just copy x1 which is already int
};

static void viFabs(MKL_INT n, const int* x, int* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = x[j] < 0 ? -x[j]: x[j];
    };
};

static void vlFabs(MKL_INT n, const long* x, long* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = x[j] < 0 ? -x[j]: x[j];
    };
};

/* Fake vsConj function just for casting purposes inside numexpr */
static void vsConj(MKL_INT n, const float* x1, float* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = x1[j];
    };
};

/* fmod not available in VML */
static void vsfmod(MKL_INT n, const float* x1, const float* x2, float* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
    dest[j] = fmodf(x1[j], x2[j]);
    };
}
static void vdfmod(MKL_INT n, const double* x1, const double* x2, double* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
    dest[j] = fmod(x1[j], x2[j]);
    };
};
// TODO: Have to add FUNC_III, FUNC_LLL signatures to functions.hpp
// static void vifmod(MKL_INT n, const int* x1, const int* x2, int* dest)
// {
//     MKL_INT j;
//     for(j=0; j < n; j++) {
//     dest[j] = fmodi(x1[j], x2[j]);
//     };
// };
// static void vlfmod(MKL_INT n, const long* x1, const long* x2, long* dest)
// {
//     MKL_INT j;
//     for(j=0; j < n; j++) {
//     dest[j] = fmodl(x1[j], x2[j]);
//     };
// };

/* no isnan, isfinite, isinf or signbit in VML */
static void vsIsfinite(MKL_INT n, const float* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isfinitef_(x1[j]);
    };
};
static void vsIsinf(MKL_INT n, const float* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isinff_(x1[j]);
    };
};
static void vsIsnan(MKL_INT n, const float* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isnanf_(x1[j]);
    };
};
static void vsSignBit(MKL_INT n, const float* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = signbitf(x1[j]);
    };
};

/* no isnan, isfinite, isinf, signbit in VML */
static void vdIsfinite(MKL_INT n, const double* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isfinited(x1[j]);
    };
};
static void vdIsinf(MKL_INT n, const double* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isinfd(x1[j]);
    };
};
static void vdIsnan(MKL_INT n, const double* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isnand(x1[j]);
    };
};
static void vdSignBit(MKL_INT n, const double* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = signbit(x1[j]);
    };
};

/* no isnan, isfinite or isinf in VML */
static void vzIsfinite(MKL_INT n, const MKL_Complex16* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isfinited(x1[j].real) && isfinited(x1[j].imag);
    };
};
static void vzIsinf(MKL_INT n, const MKL_Complex16* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isinfd(x1[j].real) || isinfd(x1[j].imag);
    };
};
static void vzIsnan(MKL_INT n, const MKL_Complex16* x1, bool* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = isnand(x1[j].real) || isnand(x1[j].imag);
    };
};

/* Fake vdConj function just for casting purposes inside numexpr */
static void vdConj(MKL_INT n, const double* x1, double* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j] = x1[j];
    };
};

/* various functions not available in VML */
static void vzExpm1(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    vzExp(n, x1, dest);
    for (j=0; j<n; j++) {
        dest[j].real -= 1.0;
    };
};

static void vzLog1p(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j].real = x1[j].real + 1;
        dest[j].imag = x1[j].imag;
    };
    vzLn(n, dest, dest);
};

static void vzLog2(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    vzLn(n, x1, dest);
    for (j=0; j<n; j++) {
        dest[j].real = dest[j].real * M_LOG2_E;
        dest[j].imag = dest[j].imag * M_LOG2_E;
    };
};

static void vzRint(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j].real = rint(x1[j].real);
        dest[j].imag = rint(x1[j].imag);
    };
};

/* Use this instead of native vzAbs in VML as it seems to work badly */
static void vzAbs_(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    for (j=0; j<n; j++) {
        dest[j].real = sqrt(x1[j].real*x1[j].real + x1[j].imag*x1[j].imag);
        dest[j].imag = 0;
    };
};

/*sign functions*/
static void vsSign(MKL_INT n, const float* x1, float* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
        dest[j] = signf(x1[j]);
    };
};
static void vdSign(MKL_INT n, const double* x1, double* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
        dest[j] = sign(x1[j]);
    };
};
static void viSign(MKL_INT n, const int* x1, int* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
        dest[j] = signi(x1[j]);
    };
};
static void vlSign(MKL_INT n, const long* x1, long* dest)
{
    MKL_INT j;
    for(j=0; j < n; j++) {
        dest[j] = signl(x1[j]);
    };
};
static void vzSign(MKL_INT n, const MKL_Complex16* x1, MKL_Complex16* dest)
{
    MKL_INT j;
    double mag;
    for(j=0; j < n; j++) {
        mag = sqrt(x1[j].real*x1[j].real + x1[j].imag*x1[j].imag);
        if (isnand(mag)) {
            dest[j].real = NAN;
            dest[j].imag = NAN;
        }
        else if (mag == 0) {
            dest[j].real = 0;
            dest[j].imag = 0;
        }
        else {
            dest[j].real = x1[j].real / mag;
            dest[j].imag = x1[j].imag / mag;
        }
    };
};
#endif
