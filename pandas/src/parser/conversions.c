
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "parser.h"

// forward declaration
double inline xstrtod(const char *p, char **q, char decimal, char sci, int skip_trailing);


inline void lowercase(char *p) {
    for ( ; *p; ++p) *p = tolower(*p);
}

inline void uppercase(char *p) {
    for ( ; *p; ++p) *p = toupper(*p);
}


/*
 *  `item` must be the nul-terminated string that is to be
 *  converted to a double.
 *
 *  To be successful, to_double() must use *all* the characters
 *  in `item`.  E.g. "1.q25" will fail.  Leading and trailing
 *  spaces are allowed.
 *
 *  `sci` is the scientific notation exponent character, usually
 *  either 'E' or 'D'.  Case is ignored.
 *
 *  `decimal` is the decimal point character, usually either
 *  '.' or ','.
 *
 */

int inline to_double(char *item, double *p_value, char sci, char decimal)
{
    char *p_end;

    *p_value = xstrtod(item, &p_end, decimal, sci, TRUE);

    return (errno == 0) && (!*p_end);
}


int inline to_complex(char *item, double *p_real, double *p_imag, char sci, char decimal)
{
    char *p_end;

    *p_real = xstrtod(item, &p_end, decimal, sci, FALSE);
    if (*p_end == '\0') {
        *p_imag = 0.0;
        return errno == 0;
    }
    if (*p_end == 'i' || *p_end == 'j') {
        *p_imag = *p_real;
        *p_real = 0.0;
        ++p_end;
    }
    else {
        if (*p_end == '+') {
            ++p_end;
        }
        *p_imag = xstrtod(p_end, &p_end, decimal, sci, FALSE);
        if (errno || ((*p_end != 'i') && (*p_end != 'j'))) {
            return FALSE;
        }
        ++p_end;
    }
    while(*p_end == ' ') {
        ++p_end;
    }
    return *p_end == '\0';
}


int inline to_longlong(char *item, long long *p_value)
{
    char *p_end;

    // Try integer conversion.  We explicitly give the base to be 10. If
    // we used 0, strtoll() would convert '012' to 10, because the leading 0 in
    // '012' signals an octal number in C.  For a general purpose reader, that
    // would be a bug, not a feature.
    *p_value = strtoll(item, &p_end, 10);

    // Allow trailing spaces.
    while (isspace(*p_end)) ++p_end;

    return (errno == 0) && (!*p_end);
}

int inline to_longlong_thousands(char *item, long long *p_value, char tsep)
{
	int i, pos, status, n = strlen(item), count = 0;
	char *tmp;
    char *p_end;

	for (i = 0; i < n; ++i)
	{
		if (*(item + i) == tsep) {
			count++;
		}
	}

	if (count == 0) {
		return to_longlong(item, p_value);
	}

	tmp = (char*) malloc((n - count + 1) * sizeof(char));
	if (tmp == NULL) {
		return 0;
	}

	pos = 0;
	for (i = 0; i < n; ++i)
	{
		if (item[i] != tsep)
			tmp[pos++] = item[i];
	}

	tmp[pos] = '\0';

	status = to_longlong(tmp, p_value);
	free(tmp);

	return status;
}

int inline to_boolean(char *item, uint8_t *val) {
	char *tmp;
	int i, status = 0;

    static const char *tstrs[2] = {"TRUE", "YES"};
    static const char *fstrs[2] = {"FALSE", "NO"};

	tmp = malloc(sizeof(char) * (strlen(item) + 1));
	strcpy(tmp, item);
	uppercase(tmp);

    for (i = 0; i < 2; ++i)
    {
        if (strcmp(tmp, tstrs[i]) == 0) {
            *val = 1;
            goto done;
        }
    }

    for (i = 0; i < 2; ++i)
    {
        if (strcmp(tmp, fstrs[i]) == 0) {
            *val = 0;
            goto done;
        }
    }

    status = -1;

done:
	free(tmp);
	return status;
}

// #define TEST

#ifdef TEST

int main(int argc, char *argv[])
{
    double x, y;
	long long xi;
    int status;
    char *s;

    //s = "0.10e-3-+5.5e2i";
    // s = "1-0j";
    // status = to_complex(s, &x, &y, 'e', '.');
	s = "123,789";
	status = to_longlong_thousands(s, &xi, ',');
    printf("s = '%s'\n", s);
    printf("status = %d\n", status);
    printf("x = %d\n", (int) xi);

    // printf("x = %lg,  y = %lg\n", x, y);

    return 0;
}
#endif

// ---------------------------------------------------------------------------
// Implementation of xstrtod

//
// strtod.c
//
// Convert string to double
//
// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// -----------------------------------------------------------------------
// Modifications by Warren Weckesser, March 2011:
// * Rename strtod() to xstrtod().
// * Added decimal and sci arguments.
// * Skip trailing spaces.
// * Commented out the other functions.
//

double inline xstrtod(const char *str, char **endptr, char decimal,
					  char sci, int skip_trailing)
{
  double number;
  int exponent;
  int negative;
  char *p = (char *) str;
  double p10;
  int n;
  int num_digits;
  int num_decimals;

  errno = 0;

  // Skip leading whitespace
  while (isspace(*p)) p++;

  // Handle optional sign
  negative = 0;
  switch (*p)
  {
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
  }

  number = 0.;
  exponent = 0;
  num_digits = 0;
  num_decimals = 0;

  // Process string of digits
  while (isdigit(*p))
  {
    number = number * 10. + (*p - '0');
    p++;
    num_digits++;
  }

  // Process decimal part
  if (*p == decimal)
  {
    p++;

    while (isdigit(*p))
    {
      number = number * 10. + (*p - '0');
      p++;
      num_digits++;
      num_decimals++;
    }

    exponent -= num_decimals;
  }

  if (num_digits == 0)
  {
    errno = ERANGE;
    return 0.0;
  }

  // Correct for sign
  if (negative) number = -number;

  // Process an exponent string
  if (toupper(*p) == toupper(sci))
  {
    // Handle optional sign
    negative = 0;
    switch (*++p)
    {
      case '-': negative = 1;   // Fall through to increment pos
      case '+': p++;
    }

    // Process string of digits
    n = 0;
    while (isdigit(*p))
    {
      n = n * 10 + (*p - '0');
      p++;
    }

    if (negative)
      exponent -= n;
    else
      exponent += n;
  }


  if (exponent < DBL_MIN_EXP  || exponent > DBL_MAX_EXP)
  {

    errno = ERANGE;
    return HUGE_VAL;
  }

  // Scale the result
  p10 = 10.;
  n = exponent;
  if (n < 0) n = -n;
  while (n)
  {
    if (n & 1)
    {
      if (exponent < 0)
        number /= p10;
      else
        number *= p10;
    }
    n >>= 1;
    p10 *= p10;
  }


  if (number == HUGE_VAL) {
	  errno = ERANGE;
  }

  if (skip_trailing) {
      // Skip trailing whitespace
      while (isspace(*p)) p++;
  }

  if (endptr) *endptr = p;


  return number;
}

/*
float strtof(const char *str, char **endptr)
{
  return (float) strtod(str, endptr);
}


long double strtold(const char *str, char **endptr)
{
  return strtod(str, endptr);
}

double atof(const char *str)
{
  return strtod(str, NULL);
}
*/

// End of xstrtod code
// ---------------------------------------------------------------------------
