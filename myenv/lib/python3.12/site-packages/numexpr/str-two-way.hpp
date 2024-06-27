/* Byte-wise substring search, using the Two-Way algorithm.
 * Copyright (C) 2008, 2010 Eric Blake
 * Permission to use, copy, modify, and distribute this software
 * is freely granted, provided that this notice is preserved.
 */


/* Before including this file, you need to include <string.h>, and define:
     RETURN_TYPE		A macro that expands to the return type.
     AVAILABLE(h, h_l, j, n_l)	A macro that returns nonzero if there are
				at least N_L bytes left starting at
				H[J].  H is 'unsigned char *', H_L, J,
				and N_L are 'size_t'; H_L is an
				lvalue.  For NUL-terminated searches,
				H_L can be modified each iteration to
				avoid having to compute the end of H
				up front.

  For case-insensitivity, you may optionally define:
     CMP_FUNC(p1, p2, l)	A macro that returns 0 iff the first L
				characters of P1 and P2 are equal.
     CANON_ELEMENT(c)		A macro that canonicalizes an element
				right after it has been fetched from
				one of the two strings.  The argument
				is an 'unsigned char'; the result must
				be an 'unsigned char' as well.

  This file undefines the macros documented above, and defines
  LONG_NEEDLE_THRESHOLD.
*/

#include <limits.h>

/*
  Python 2.7 (the only Python 2.x version supported as of now and until 2020)
  is built on windows with Visual Studio 2008 C compiler. That dictates that
  the compiler which must be used by authors of third party Python modules.
  See https://mail.python.org/pipermail/distutils-sig/2014-September/024885.html

  Unfortunately this version of Visual Studio doesn't claim to be C99 compatible
  and in particular it lacks the stdint.h header. So we have to replace it with
  a public domain version.

  Visual Studio 2010 and later have stdint.h.
*/

#ifdef _MSC_VER
  #if _MSC_VER <= 1500
    #include "win32/stdint.h"
  #endif
#else
  #include <stdint.h>
#endif

/* We use the Two-Way string matching algorithm, which guarantees
   linear complexity with constant space.  Additionally, for long
   needles, we also use a bad character shift table similar to the
   Boyer-Moore algorithm to achieve improved (potentially sub-linear)
   performance.

   See http://www-igm.univ-mlv.fr/~lecroq/string/node26.html#SECTION00260
   and http://en.wikipedia.org/wiki/Boyer-Moore_string_search_algorithm
*/

/* Point at which computing a bad-byte shift table is likely to be
   worthwhile.  Small needles should not compute a table, since it
   adds (1 << CHAR_BIT) + NEEDLE_LEN computations of preparation for a
   speedup no greater than a factor of NEEDLE_LEN.  The larger the
   needle, the better the potential performance gain.  On the other
   hand, on non-POSIX systems with CHAR_BIT larger than eight, the
   memory required for the table is prohibitive.  */
#if CHAR_BIT < 10
# define LONG_NEEDLE_THRESHOLD 32U
#else
# define LONG_NEEDLE_THRESHOLD SIZE_MAX
#endif

#define MAX(a, b) ((a < b) ? (b) : (a))

#ifndef CANON_ELEMENT
# define CANON_ELEMENT(c) c
#endif
#ifndef CMP_FUNC
# define CMP_FUNC memcmp
#endif

/* Perform a critical factorization of NEEDLE, of length NEEDLE_LEN.
   Return the index of the first byte in the right half, and set
   *PERIOD to the global period of the right half.

   The global period of a string is the smallest index (possibly its
   length) at which all remaining bytes in the string are repetitions
   of the prefix (the last repetition may be a subset of the prefix).

   When NEEDLE is factored into two halves, a local period is the
   length of the smallest word that shares a suffix with the left half
   and shares a prefix with the right half.  All factorizations of a
   non-empty NEEDLE have a local period of at least 1 and no greater
   than NEEDLE_LEN.

   A critical factorization has the property that the local period
   equals the global period.  All strings have at least one critical
   factorization with the left half smaller than the global period.

   Given an ordered alphabet, a critical factorization can be computed
   in linear time, with 2 * NEEDLE_LEN comparisons, by computing the
   larger of two ordered maximal suffixes.  The ordered maximal
   suffixes are determined by lexicographic comparison of
   periodicity.  */
static size_t
critical_factorization (const unsigned char *needle, size_t needle_len,
			size_t *period)
{
  /* Index of last byte of left half, or SIZE_MAX.  */
  size_t max_suffix, max_suffix_rev;
  size_t j; /* Index into NEEDLE for current candidate suffix.  */
  size_t k; /* Offset into current period.  */
  size_t p; /* Intermediate period.  */
  unsigned char a, b; /* Current comparison bytes.  */

  /* Invariants:
     0 <= j < NEEDLE_LEN - 1
     -1 <= max_suffix{,_rev} < j (treating SIZE_MAX as if it were signed)
     min(max_suffix, max_suffix_rev) < global period of NEEDLE
     1 <= p <= global period of NEEDLE
     p == global period of the substring NEEDLE[max_suffix{,_rev}+1...j]
     1 <= k <= p
  */

  /* Perform lexicographic search.  */
  max_suffix = SIZE_MAX;
  j = 0;
  k = p = 1;
  while (j + k < needle_len)
    {
      a = CANON_ELEMENT (needle[j + k]);
      b = CANON_ELEMENT (needle[(size_t)(max_suffix + k)]);
      if (a < b)
	{
	  /* Suffix is smaller, period is entire prefix so far.  */
	  j += k;
	  k = 1;
	  p = j - max_suffix;
	}
      else if (a == b)
	{
	  /* Advance through repetition of the current period.  */
	  if (k != p)
	    ++k;
	  else
	    {
	      j += p;
	      k = 1;
	    }
	}
      else /* b < a */
	{
	  /* Suffix is larger, start over from current location.  */
	  max_suffix = j++;
	  k = p = 1;
	}
    }
  *period = p;

  /* Perform reverse lexicographic search.  */
  max_suffix_rev = SIZE_MAX;
  j = 0;
  k = p = 1;
  while (j + k < needle_len)
    {
      a = CANON_ELEMENT (needle[j + k]);
      b = CANON_ELEMENT (needle[max_suffix_rev + k]);
      if (b < a)
	{
	  /* Suffix is smaller, period is entire prefix so far.  */
	  j += k;
	  k = 1;
	  p = j - max_suffix_rev;
	}
      else if (a == b)
	{
	  /* Advance through repetition of the current period.  */
	  if (k != p)
	    ++k;
	  else
	    {
	      j += p;
	      k = 1;
	    }
	}
      else /* a < b */
	{
	  /* Suffix is larger, start over from current location.  */
	  max_suffix_rev = j++;
	  k = p = 1;
	}
    }

  /* Choose the longer suffix.  Return the first byte of the right
     half, rather than the last byte of the left half.  */
  if (max_suffix_rev + 1 < max_suffix + 1)
    return max_suffix + 1;
  *period = p;
  return max_suffix_rev + 1;
}

/* Return the first location of non-empty NEEDLE within HAYSTACK, or
   NULL.  HAYSTACK_LEN is the minimum known length of HAYSTACK.  This
   method is optimized for NEEDLE_LEN < LONG_NEEDLE_THRESHOLD.
   Performance is guaranteed to be linear, with an initialization cost
   of 2 * NEEDLE_LEN comparisons.

   If AVAILABLE does not modify HAYSTACK_LEN (as in memmem), then at
   most 2 * HAYSTACK_LEN - NEEDLE_LEN comparisons occur in searching.
   If AVAILABLE modifies HAYSTACK_LEN (as in strstr), then at most 3 *
   HAYSTACK_LEN - NEEDLE_LEN comparisons occur in searching.  */
static RETURN_TYPE
two_way_short_needle (const unsigned char *haystack, size_t haystack_len,
		      const unsigned char *needle, size_t needle_len)
{
  size_t i; /* Index into current byte of NEEDLE.  */
  size_t j; /* Index into current window of HAYSTACK.  */
  size_t period; /* The period of the right half of needle.  */
  size_t suffix; /* The index of the right half of needle.  */

  /* Factor the needle into two halves, such that the left half is
     smaller than the global period, and the right half is
     periodic (with a period as large as NEEDLE_LEN - suffix).  */
  suffix = critical_factorization (needle, needle_len, &period);

  /* Perform the search.  Each iteration compares the right half
     first.  */
  if (CMP_FUNC (needle, needle + period, suffix) == 0)
    {
      /* Entire needle is periodic; a mismatch can only advance by the
	 period, so use memory to avoid rescanning known occurrences
	 of the period.  */
      size_t memory = 0;
      j = 0;
      while (AVAILABLE (haystack, haystack_len, j, needle_len))
	{
	  /* Scan for matches in right half.  */
	  i = MAX (suffix, memory);
	  while (i < needle_len && (CANON_ELEMENT (needle[i])
				    == CANON_ELEMENT (haystack[i + j])))
	    ++i;
	  if (needle_len <= i)
	    {
	      /* Scan for matches in left half.  */
	      i = suffix - 1;
	      while (memory < i + 1 && (CANON_ELEMENT (needle[i])
					== CANON_ELEMENT (haystack[i + j])))
		--i;
	      if (i + 1 < memory + 1)
		return (RETURN_TYPE) (haystack + j);
	      /* No match, so remember how many repetitions of period
		 on the right half were scanned.  */
	      j += period;
	      memory = needle_len - period;
	    }
	  else
	    {
	      j += i - suffix + 1;
	      memory = 0;
	    }
	}
    }
  else
    {
      /* The two halves of needle are distinct; no extra memory is
	 required, and any mismatch results in a maximal shift.  */
      period = MAX (suffix, needle_len - suffix) + 1;
      j = 0;
      while (AVAILABLE (haystack, haystack_len, j, needle_len))
	{
	  /* Scan for matches in right half.  */
	  i = suffix;
	  while (i < needle_len && (CANON_ELEMENT (needle[i])
				    == CANON_ELEMENT (haystack[i + j])))
	    ++i;
	  if (needle_len <= i)
	    {
	      /* Scan for matches in left half.  */
	      i = suffix - 1;
	      while (i != SIZE_MAX && (CANON_ELEMENT (needle[i])
				       == CANON_ELEMENT (haystack[i + j])))
		--i;
	      if (i == SIZE_MAX)
		return (RETURN_TYPE) (haystack + j);
	      j += period;
	    }
	  else
	    j += i - suffix + 1;
	}
    }
  return NULL;
}

/* Return the first location of non-empty NEEDLE within HAYSTACK, or
   NULL.  HAYSTACK_LEN is the minimum known length of HAYSTACK.  This
   method is optimized for LONG_NEEDLE_THRESHOLD <= NEEDLE_LEN.
   Performance is guaranteed to be linear, with an initialization cost
   of 3 * NEEDLE_LEN + (1 << CHAR_BIT) operations.

   If AVAILABLE does not modify HAYSTACK_LEN (as in memmem), then at
   most 2 * HAYSTACK_LEN - NEEDLE_LEN comparisons occur in searching,
   and sublinear performance O(HAYSTACK_LEN / NEEDLE_LEN) is possible.
   If AVAILABLE modifies HAYSTACK_LEN (as in strstr), then at most 3 *
   HAYSTACK_LEN - NEEDLE_LEN comparisons occur in searching, and
   sublinear performance is not possible.  */
static RETURN_TYPE
two_way_long_needle (const unsigned char *haystack, size_t haystack_len,
		     const unsigned char *needle, size_t needle_len)
{
  size_t i; /* Index into current byte of NEEDLE.  */
  size_t j; /* Index into current window of HAYSTACK.  */
  size_t period; /* The period of the right half of needle.  */
  size_t suffix; /* The index of the right half of needle.  */
  size_t shift_table[1U << CHAR_BIT]; /* See below.  */

  /* Factor the needle into two halves, such that the left half is
     smaller than the global period, and the right half is
     periodic (with a period as large as NEEDLE_LEN - suffix).  */
  suffix = critical_factorization (needle, needle_len, &period);

  /* Populate shift_table.  For each possible byte value c,
     shift_table[c] is the distance from the last occurrence of c to
     the end of NEEDLE, or NEEDLE_LEN if c is absent from the NEEDLE.
     shift_table[NEEDLE[NEEDLE_LEN - 1]] contains the only 0.  */
  for (i = 0; i < 1U << CHAR_BIT; i++)
    shift_table[i] = needle_len;
  for (i = 0; i < needle_len; i++)
    shift_table[CANON_ELEMENT (needle[i])] = needle_len - i - 1;

  /* Perform the search.  Each iteration compares the right half
     first.  */
  if (CMP_FUNC (needle, needle + period, suffix) == 0)
    {
      /* Entire needle is periodic; a mismatch can only advance by the
	 period, so use memory to avoid rescanning known occurrences
	 of the period.  */
      size_t memory = 0;
      size_t shift;
      j = 0;
      while (AVAILABLE (haystack, haystack_len, j, needle_len))
	{
	  /* Check the last byte first; if it does not match, then
	     shift to the next possible match location.  */
	  shift = shift_table[CANON_ELEMENT (haystack[j + needle_len - 1])];
	  if (0 < shift)
	    {
	      if (memory && shift < period)
		{
		  /* Since needle is periodic, but the last period has
		     a byte out of place, there can be no match until
		     after the mismatch.  */
		  shift = needle_len - period;
		}
	      memory = 0;
	      j += shift;
	      continue;
	    }
	  /* Scan for matches in right half.  The last byte has
	     already been matched, by virtue of the shift table.  */
	  i = MAX (suffix, memory);
	  while (i < needle_len - 1 && (CANON_ELEMENT (needle[i])
					== CANON_ELEMENT (haystack[i + j])))
	    ++i;
	  if (needle_len - 1 <= i)
	    {
	      /* Scan for matches in left half.  */
	      i = suffix - 1;
	      while (memory < i + 1 && (CANON_ELEMENT (needle[i])
					== CANON_ELEMENT (haystack[i + j])))
		--i;
	      if (i + 1 < memory + 1)
		return (RETURN_TYPE) (haystack + j);
	      /* No match, so remember how many repetitions of period
		 on the right half were scanned.  */
	      j += period;
	      memory = needle_len - period;
	    }
	  else
	    {
	      j += i - suffix + 1;
	      memory = 0;
	    }
	}
    }
  else
    {
      /* The two halves of needle are distinct; no extra memory is
	 required, and any mismatch results in a maximal shift.  */
      size_t shift;
      period = MAX (suffix, needle_len - suffix) + 1;
      j = 0;
      while (AVAILABLE (haystack, haystack_len, j, needle_len))
	{
	  /* Check the last byte first; if it does not match, then
	     shift to the next possible match location.  */
	  shift = shift_table[CANON_ELEMENT (haystack[j + needle_len - 1])];
	  if (0 < shift)
	    {
	      j += shift;
	      continue;
	    }
	  /* Scan for matches in right half.  The last byte has
	     already been matched, by virtue of the shift table.  */
	  i = suffix;
	  while (i < needle_len - 1 && (CANON_ELEMENT (needle[i])
					== CANON_ELEMENT (haystack[i + j])))
	    ++i;
	  if (needle_len - 1 <= i)
	    {
	      /* Scan for matches in left half.  */
	      i = suffix - 1;
	      while (i != SIZE_MAX && (CANON_ELEMENT (needle[i])
				       == CANON_ELEMENT (haystack[i + j])))
		--i;
	      if (i == SIZE_MAX)
		return (RETURN_TYPE) (haystack + j);
	      j += period;
	    }
	  else
	    j += i - suffix + 1;
	}
    }
  return NULL;
}

#undef AVAILABLE
#undef CANON_ELEMENT
#undef CMP_FUNC
#undef MAX
#undef RETURN_TYPE
