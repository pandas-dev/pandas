#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

#include "libbase64.h"
#include "../../tables/tables.h"
#include "../../codecs.h"
#include "config.h"
#include "../../env.h"

#if HAVE_AVX512
#if defined(__clang__)
#pragma clang attribute push (__attribute__((target("avx512vbmi"))), apply_to=function)
#else
#pragma GCC target("avx512vbmi")
#endif
#include <immintrin.h>

#include "../avx2/dec_reshuffle.c"
#include "../avx2/dec_loop.c"
#include "enc_reshuffle_translate.c"
#include "enc_loop.c"

#endif	// HAVE_AVX512

void
base64_stream_encode_avx512 BASE64_ENC_PARAMS
{
#if HAVE_AVX512
	#include "../generic/enc_head.c"
	enc_loop_avx512(&s, &slen, &o, &olen);
	#include "../generic/enc_tail.c"
#else
	base64_enc_stub(state, src, srclen, out, outlen);
#endif
}

// Reuse AVX2 decoding. Not supporting AVX512 at present
int
base64_stream_decode_avx512 BASE64_DEC_PARAMS
{
#if HAVE_AVX512
	#include "../generic/dec_head.c"
	dec_loop_avx2(&s, &slen, &o, &olen);
	#include "../generic/dec_tail.c"
#if defined(__clang__)
	#pragma clang attribute pop
#endif
#else
	return base64_dec_stub(state, src, srclen, out, outlen);
#endif
}
