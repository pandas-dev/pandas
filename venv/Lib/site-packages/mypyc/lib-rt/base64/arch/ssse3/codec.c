#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#include "libbase64.h"
#include "../../tables/tables.h"
#include "../../codecs.h"
#include "config.h"
#include "../../env.h"

#if HAVE_SSSE3
#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("ssse3"))), apply_to=function)
#else
#pragma GCC target("ssse3")
#endif
#include <tmmintrin.h>

// Only enable inline assembly on supported compilers and on 64-bit CPUs.
// 32-bit CPUs with SSSE3 support, such as low-end Atoms, only have eight XMM
// registers, which is not enough to run the inline assembly.
#ifndef BASE64_SSSE3_USE_ASM
# if (defined(__GNUC__) || defined(__clang__)) && BASE64_WORDSIZE == 64
#  define BASE64_SSSE3_USE_ASM 1
# else
#  define BASE64_SSSE3_USE_ASM 0
# endif
#endif

#include "dec_reshuffle.c"
#include "dec_loop.c"

#if BASE64_SSSE3_USE_ASM
# include "enc_loop_asm.c"
#else
# include "enc_reshuffle.c"
# include "enc_translate.c"
# include "enc_loop.c"
#endif

#endif	// HAVE_SSSE3

void
base64_stream_encode_ssse3 BASE64_ENC_PARAMS
{
#if HAVE_SSSE3
	#include "../generic/enc_head.c"
	enc_loop_ssse3(&s, &slen, &o, &olen);
	#include "../generic/enc_tail.c"
#else
	base64_enc_stub(state, src, srclen, out, outlen);
#endif
}

int
base64_stream_decode_ssse3 BASE64_DEC_PARAMS
{
#if HAVE_SSSE3
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
