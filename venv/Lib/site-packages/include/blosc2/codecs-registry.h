/*********************************************************************
  Blosc - Blocked Shuffling and Compression Library

  Copyright (c) 2021  Blosc Development Team <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

#ifndef BLOSC_BLOSC2_CODECS_REGISTRY_H
#define BLOSC_BLOSC2_CODECS_REGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
    BLOSC_CODEC_NDLZ = 32,
    //!< Simple Lempel-Ziv compressor for NDim data. Experimental, mainly for teaching purposes.
    BLOSC_CODEC_ZFP_FIXED_ACCURACY = 33,
    //!< ZFP compressor for fixed accuracy mode. The desired accuracy is set in `compcode_meta`.
    //!< See https://github.com/Blosc/c-blosc2/blob/main/plugins/codecs/zfp/README.md
    BLOSC_CODEC_ZFP_FIXED_PRECISION = 34,
    //!< ZFP compressor for fixed precision. The desired precision is set in `compcode_meta`.
    //!< See https://github.com/Blosc/c-blosc2/blob/main/plugins/codecs/zfp/README.md
    BLOSC_CODEC_ZFP_FIXED_RATE = 35,
    //!< ZFP compressor for fixed precision. The desired rate is set in `compcode_meta`.
    //!< See https://github.com/Blosc/c-blosc2/blob/main/plugins/codecs/zfp/README.md
    BLOSC_CODEC_OPENHTJ2K = 36,
    //!< OpenHTJ2K compressor for JPEG 2000 HT.
    //!< See https://github.com/Blosc/blosc2_openhtj2k
    BLOSC_CODEC_GROK = 37,
    //!< Grok compressor for JPEG 2000.
    //!< See https://github.com/Blosc/blosc2_grok
};

void register_codecs(void);

// For dynamically loaded codecs
typedef struct {
    char *encoder;
    char *decoder;
} codec_info;

#ifdef __cplusplus
}
#endif

#endif /* BLOSC_BLOSC2_CODECS_REGISTRY_H */
