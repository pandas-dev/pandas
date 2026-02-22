#ifndef NUMEXPR_COMPLEX_FUNCTIONS_HPP
#define NUMEXPR_COMPLEX_FUNCTIONS_HPP

/*********************************************************************
  Numexpr - Fast numerical array expression evaluator for NumPy.

      License: MIT
      Author:  See AUTHORS.txt

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

// Replace npy_cdouble with std::complex<double>
#include <math.h> // NAN
#include <complex>

/* constants */
static std::complex<double> nc_1(1., 0.);
static std::complex<double> nc_half(0.5, 0.);
static std::complex<double> nc_i(0., 1.);
static std::complex<double> nc_i2(0., 0.5);
/*
static std::complex<double> nc_mi = {0., -1.};
static std::complex<double> nc_pi2 = {M_PI/2., 0.};
*/

/* *************************** WARNING *****************************
Due to the way Numexpr places the results of operations, the *x and *r
pointers do point to the same address (apparently this doesn't happen
in NumPy).  So, measures should be taken so as to not to reuse *x
after the first *r has been overwritten.
*********************************************************************
*/

static void
nc_assign(std::complex<double> *x, std::complex<double> *r)
{
  r->real(x->real());
  r->imag(x->imag());
  return;
}

static void
nc_sum(std::complex<double> *a, std::complex<double> *b, std::complex<double> *r)
{
    r->real(a->real() + b->real());
    r->imag(a->imag() + b->imag());
    return;
}

static void
nc_diff(std::complex<double> *a, std::complex<double> *b, std::complex<double> *r)
{
    r->real(a->real() - b->real());
    r->imag(a->imag() - b->imag());
    return;
}

static void
nc_neg(std::complex<double> *a, std::complex<double> *r)
{
    r->real(-a->real());
    r->imag(-a->imag());
    return;
}

static void
nc_conj(std::complex<double> *a, std::complex<double> *r)
{
    r->real(a->real());
    r->imag(-a->imag());
    return;
}

// Needed for allowing the internal casting in numexpr machinery for
// conjugate operations
inline float fconjf(float x)
{
    return x;
}

// Needed for allowing the internal casting in numexpr machinery for
// conjugate operations
inline double fconj(double x)
{
    return x;
}

static void
nc_prod(std::complex<double> *a, std::complex<double> *b, std::complex<double> *r)
{
    double ar=a->real(), br=b->real(), ai=a->imag(), bi=b->imag();
    r->real(ar*br - ai*bi);
    r->imag(ar*bi + ai*br);
    return;
}

static void
nc_quot(std::complex<double> *a, std::complex<double> *b, std::complex<double> *r)
{
    double ar=a->real(), br=b->real(), ai=a->imag(), bi=b->imag();
    double d = br*br + bi*bi;
    r->real((ar*br + ai*bi)/d);
    r->imag((ai*br - ar*bi)/d);
    return;
}

static void
nc_sqrt(std::complex<double> *x, std::complex<double> *r)
{
    double s,d;
    if (x->real() == 0. && x->imag() == 0.)
        *r = *x;
    else {
        s = sqrt((fabs(x->real()) + hypot(x->real(),x->imag()))/2);
        d = x->imag()/(2*s);
        if (x->real() > 0.) {
            r->real(s);
            r->imag(d);
        }
        else if (x->imag() >= 0.) {
            r->real(d);
            r->imag(s);
        }
        else {
            r->real(-d);
            r->imag(-s);
        }
    }
    return;
}

static void
nc_log(std::complex<double> *x, std::complex<double> *r)
{
    double l = hypot(x->real(),x->imag());
    r->imag(atan2(x->imag(), x->real()));
    r->real(log(l));
    return;
}

static void
nc_log1p(std::complex<double> *x, std::complex<double> *r)
{
    double l = hypot(x->real() + 1.0,x->imag());
    r->imag(atan2(x->imag(), x->real() + 1.0));
    r->real(log(l));
    return;
}

static void
nc_exp(std::complex<double> *x, std::complex<double> *r)
{
    double a = exp(x->real());
    r->real(a*cos(x->imag()));
    r->imag(a*sin(x->imag()));
    return;
}

static void
nc_expm1(std::complex<double> *x, std::complex<double> *r)
{
    double a = sin(x->imag() / 2);
    double b = exp(x->real());
    r->real(expm1(x->real()) * cos(x->imag()) - 2 * a * a);
    r->imag(b * sin(x->imag()));
    return;
}

static void
nc_pow(std::complex<double> *a, std::complex<double> *b, std::complex<double> *r)
{
    npy_intp n;
    double ar=a->real(), br=b->real(), ai=a->imag(), bi=b->imag();

    if (br == 0. && bi == 0.) {
        r->real(1.);
        r->imag(0.);
        return;
    }
    if (ar == 0. && ai == 0.) {
        r->real(0.);
        r->imag(0.);
        return;
    }
    if (bi == 0 && (n=(npy_intp)br) == br) {
        if (n > -100 && n < 100) {
        std::complex<double> p, aa;
        npy_intp mask = 1;
        if (n < 0) n = -n;
        aa = nc_1;
        p.real(ar); p.imag(ai);
        while (1) {
            if (n & mask)
                nc_prod(&aa,&p,&aa);
            mask <<= 1;
            if (n < mask || mask <= 0) break;
            nc_prod(&p,&p,&p);
        }
        r->real(aa.real()); r->imag(aa.imag());
        if (br < 0) nc_quot(&nc_1, r, r);
        return;
        }
    }
    /* complexobject.c uses an inline version of this formula
       investigate whether this had better performance or accuracy */
    nc_log(a, r);
    nc_prod(r, b, r);
    nc_exp(r, r);
    return;
}


static void
nc_prodi(std::complex<double> *x, std::complex<double> *r)
{
    double xr = x->real();
    r->real(-x->imag());
    r->imag(xr);
    return;
}


static void
nc_acos(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> a, *pa=&a;

    nc_assign(x, pa);
    nc_prod(x,x,r);
    nc_diff(&nc_1, r, r);
    nc_sqrt(r, r);
    nc_prodi(r, r);
    nc_sum(pa, r, r);
    nc_log(r, r);
    nc_prodi(r, r);
    nc_neg(r, r);
    return;
    /* return nc_neg(nc_prodi(nc_log(nc_sum(x,nc_prod(nc_i,
       nc_sqrt(nc_diff(nc_1,nc_prod(x,x))))))));
    */
}

static void
nc_acosh(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> t, a, *pa=&a;

    nc_assign(x, pa);
    nc_sum(x, &nc_1, &t);
    nc_sqrt(&t, &t);
    nc_diff(x, &nc_1, r);
    nc_sqrt(r, r);
    nc_prod(&t, r, r);
    nc_sum(pa, r, r);
    nc_log(r, r);
    return;
    /*
      return nc_log(nc_sum(x,
      nc_prod(nc_sqrt(nc_sum(x,nc_1)), nc_sqrt(nc_diff(x,nc_1)))));
    */
}

static void
nc_asin(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> a, *pa=&a;
    nc_prodi(x, pa);
    nc_prod(x, x, r);
    nc_diff(&nc_1, r, r);
    nc_sqrt(r, r);
    nc_sum(pa, r, r);
    nc_log(r, r);
    nc_prodi(r, r);
    nc_neg(r, r);
    return;
    /*
      return nc_neg(nc_prodi(nc_log(nc_sum(nc_prod(nc_i,x),
      nc_sqrt(nc_diff(nc_1,nc_prod(x,x)))))));
    */
}


static void
nc_asinh(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> a, *pa=&a;
    nc_assign(x, pa);
    nc_prod(x, x, r);
    nc_sum(&nc_1, r, r);
    nc_sqrt(r, r);
    nc_sum(r, pa, r);
    nc_log(r, r);
    return;
    /*
      return nc_log(nc_sum(nc_sqrt(nc_sum(nc_1,nc_prod(x,x))),x));
    */
}

static void
nc_atan(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> a, *pa=&a;
    nc_diff(&nc_i, x, pa);
    nc_sum(&nc_i, x, r);
    nc_quot(r, pa, r);
    nc_log(r,r);
    nc_prod(&nc_i2, r, r);
    return;
    /*
      return nc_prod(nc_i2,nc_log(nc_quot(nc_sum(nc_i,x),nc_diff(nc_i,x))));
    */
}

static void
nc_atanh(std::complex<double> *x, std::complex<double> *r)
{
    std::complex<double> a, b, *pa=&a, *pb=&b;
    nc_assign(x, pa);
    nc_diff(&nc_1, pa, r);
    nc_sum(&nc_1, pa, pb);
    nc_quot(pb, r, r);
    nc_log(r, r);
    nc_prod(&nc_half, r, r);
    return;
    /*
      return nc_prod(nc_half,nc_log(nc_quot(nc_sum(nc_1,x),nc_diff(nc_1,x))));
    */
}

static void
nc_cos(std::complex<double> *x, std::complex<double> *r)
{
    double xr=x->real(), xi=x->imag();
    r->real(cos(xr)*cosh(xi));
    r->imag(-sin(xr)*sinh(xi));
    return;
}

static void
nc_cosh(std::complex<double> *x, std::complex<double> *r)
{
    double xr=x->real(), xi=x->imag();
    r->real(cos(xi)*cosh(xr));
    r->imag(sin(xi)*sinh(xr));
    return;
}


#define M_LOG10_E 0.434294481903251827651128918916605082294397
#define M_LOG2_E  1.44269504088896340735992468100189213742664


static void
nc_log10(std::complex<double> *x, std::complex<double> *r)
{
    nc_log(x, r);
    r->real(r->real() * M_LOG10_E);
    r->imag(r->imag() * M_LOG10_E);
    return;
}

static void
nc_log2(std::complex<double> *x, std::complex<double> *r)
{
    nc_log(x, r);
    r->real(r->real() * M_LOG2_E);
    r->imag(r->imag() * M_LOG2_E);
    return;
}

static void
nc_sin(std::complex<double> *x, std::complex<double> *r)
{
    double xr=x->real(), xi=x->imag();
    r->real(sin(xr)*cosh(xi));
    r->imag(cos(xr)*sinh(xi));
    return;
}

static void
nc_sinh(std::complex<double> *x, std::complex<double> *r)
{
    double xr=x->real(), xi=x->imag();
    r->real(cos(xi)*sinh(xr));
    r->imag(sin(xi)*cosh(xr));
    return;
}

static void
nc_tan(std::complex<double> *x, std::complex<double> *r)
{
    double xr = x->real();
    double xi = x->imag();
    double imag_part;

    double denom = cos(2*xr) + cosh(2*xi);
    // handle overflows
    if (xi > 20) {
        imag_part = 1.0 / (1.0 + exp(-4*xi));
    } else if (xi < -20) {
        imag_part = -1.0 / (1.0 + exp(4*xi));
    } else {
        imag_part = sinh(2*xi) / denom;
    }
    double real_part = sin(2*xr) / denom;

    r->real(real_part);
    r->imag(imag_part);
    return;
}

static void
nc_tanh(std::complex<double> *x, std::complex<double> *r)
{
    double xr = x->real();
    double xi = x->imag();
    double real_part;
    double denom = cosh(2*xr) + cos(2*xi);
    // handle overflows
    if (xr > 20) {
        real_part = 1.0 / (1.0 + exp(-4*xr));
    } else if (xr < -20) {
        real_part = -1.0 / (1.0 + exp(4*xr));
    } else {
        real_part = sinh(2*xr) / denom;
    }
    double imag_part = sin(2*xi) / denom;

    r->real(real_part);
    r->imag(imag_part);
    return;
}

static void
nc_abs(std::complex<double> *x, std::complex<double> *r)
{
    r->real(sqrt(x->real()*x->real() + x->imag()*x->imag()));
    r->imag(0);
}

static void
nc_rint(std::complex<double> *x, std::complex<double> *r)
{
    r->real(rint(x->real()));
    r->imag(rint(x->imag()));
}

static bool
nc_isinf(std::complex<double> *x)
{
    double xr=x->real(), xi=x->imag();
    bool bi,br;
    bi = isinfd(xi);
    br = isinfd(xr);
    return bi || br;
}

static bool
nc_isnan(std::complex<double> *x)
{
    double xr=x->real(), xi=x->imag();
    bool bi,br;
    bi = isnand(xi);
    br = isnand(xr);
    return bi || br;
}

static bool
nc_isfinite(std::complex<double> *x)
{
    double xr=x->real(), xi=x->imag();
    bool bi,br;
    bi = isfinited(xi);
    br = isfinited(xr);
    return bi && br;
}

static void
nc_sign(std::complex<double> *x, std::complex<double> *r)
{
    if (nc_isnan(x)){
        r->real(NAN);
        r->imag(NAN);
    }
    std::complex<double> mag;
    nc_abs(x, &mag);
    if (mag.real() == 0){
        r->real(0);
        r->imag(0);
    }
    else{
        r->real(x->real()/mag.real());
        r->imag(x->imag()/mag.real());
    }
}

#endif // NUMEXPR_COMPLEX_FUNCTIONS_HPP
