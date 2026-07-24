#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#include "libbase64.h"
#include "../../tables/tables.h"
#include "../../codecs.h"
#include "config.h"
#include "../../env.h"

#if HAVE_SSE42
#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("sse4.2"))), apply_to=function)
#else
#pragma GCC target("sse4.2")
#endif
#include <nmmintrin.h>

// Only enable inline assembly on supported compilers and on 64-bit CPUs.
#ifndef BASE64_SSE42_USE_ASM
# if (defined(__GNUC__) || defined(__clang__)) && BASE64_WORDSIZE == 64
#  define BASE64_SSE42_USE_ASM 1
# else
#  define BASE64_SSE42_USE_ASM 0
# endif
#endif

#include "../ssse3/dec_reshuffle.c"
#include "../ssse3/dec_loop.c"

#if BASE64_SSE42_USE_ASM
# include "../ssse3/enc_loop_asm.c"
#else
# include "../ssse3/enc_translate.c"
# include "../ssse3/enc_reshuffle.c"
# include "../ssse3/enc_loop.c"
#endif

#endif	// HAVE_SSE42

void
base64_stream_encode_sse42 BASE64_ENC_PARAMS
{
#if HAVE_SSE42
	#include "../generic/enc_head.c"
	enc_loop_ssse3(&s, &slen, &o, &olen);
	#include "../generic/enc_tail.c"
#else
	base64_enc_stub(state, src, srclen, out, outlen);
#endif
}

int
base64_stream_decode_sse42 BASE64_DEC_PARAMS
{
#if HAVE_SSE42
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
