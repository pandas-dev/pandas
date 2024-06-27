/*********************************************************************
  Blosc - Blocked Shuffling and Compression Library

  Copyright (C) 2021  The Blosc Developers <blosc@blosc.org>
  https://blosc.org
  License: BSD 3-Clause (see LICENSE.txt)

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

#ifndef BLOSC_BLOSC2_TUNERS_REGISTRY_H
#define BLOSC_BLOSC2_TUNERS_REGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
    BLOSC_BTUNE = 32,
};

void register_tuners(void);

// For dynamically loaded tuners
typedef struct {
    char *init;
    char *next_blocksize;
    char *next_cparams;
    char *update;
    char *free;
} tuner_info;

#ifdef __cplusplus
}
#endif

#endif /* BLOSC_BLOSC2_TUNERS_REGISTRY_H */
