#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#include "libbase64.h"
#include "../../tables/tables.h"
#include "../../codecs.h"
#include "config.h"
#include "../../env.h"

#if HAVE_AVX
#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("avx"))), apply_to=function)
#else
#pragma GCC target("avx")
#endif
#include <immintrin.h>

// Only enable inline assembly on supported compilers and on 64-bit CPUs.
#ifndef BASE64_AVX_USE_ASM
# if (defined(__GNUC__) || defined(__clang__)) && BASE64_WORDSIZE == 64
#  define BASE64_AVX_USE_ASM 1
# else
#  define BASE64_AVX_USE_ASM 0
# endif
#endif

#include "../ssse3/dec_reshuffle.c"
#include "../ssse3/dec_loop.c"

#if BASE64_AVX_USE_ASM
# include "./enc_loop_asm.c"
#else
# include "../ssse3/enc_translate.c"
# include "../ssse3/enc_reshuffle.c"
# include "../ssse3/enc_loop.c"
#endif

#endif	// HAVE_AVX

void
base64_stream_encode_avx BASE64_ENC_PARAMS
{
#if HAVE_AVX
	#include "../generic/enc_head.c"

	// For supported compilers, use a hand-optimized inline assembly
	// encoder. Otherwise fall back on the SSSE3 encoder, but compiled with
	// AVX flags to generate better optimized AVX code.

#if BASE64_AVX_USE_ASM
	enc_loop_avx(&s, &slen, &o, &olen);
#else
	enc_loop_ssse3(&s, &slen, &o, &olen);
#endif

	#include "../generic/enc_tail.c"
#else
	base64_enc_stub(state, src, srclen, out, outlen);
#endif
}

int
base64_stream_decode_avx BASE64_DEC_PARAMS
{
#if HAVE_AVX
	#include "../generic/dec_head.c"
	dec_loop_ssse3(&s, &slen, &o, &olen);
	#include "../generic/dec_tail.c"
#if defined(__clang__)
	#pragma clang attribute pop
#endif
#else
	return base64_dec_stub(state, src, srclen, out, outlen);
#endif
}
