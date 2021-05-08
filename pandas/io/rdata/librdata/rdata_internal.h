/*
Copyright (c) 2020 Evan Miller
*/

//
//  rdata_internal.h
//

#ifndef PANDAS_IO_RDATA_LIBRDATA_RDATA_INTERNAL_H_
#define PANDAS_IO_RDATA_LIBRDATA_RDATA_INTERNAL_H_

#include "rdata_bits.h"

#pragma pack(push, 1)

typedef struct rdata_v2_header_s {
    char       header[2];
    uint32_t   format_version;
    uint32_t   writer_version;
    uint32_t   reader_version;
} rdata_v2_header_t;

typedef struct rdata_sexptype_header_s {
    unsigned int   type:8;
    unsigned int   object:1;
    unsigned int   attributes:1;
    unsigned int   tag:1;
    unsigned int   unused:1;
    unsigned int   gp:16;
    unsigned int   padding:4;
} rdata_sexptype_header_t;

typedef struct rdata_sexptype_info_s {
    rdata_sexptype_header_t  header;
    int32_t                  attributes;
    int32_t                  tag;
    int32_t                  ref;
} rdata_sexptype_info_t;

#pragma pack(pop)

#define RDATA_SEXPTYPE_NIL                 0
#define RDATA_SEXPTYPE_SYMBOL              1
#define RDATA_SEXPTYPE_PAIRLIST            2
#define RDATA_SEXPTYPE_CLOSURE             3
#define RDATA_SEXPTYPE_ENVIRONMENT         4
#define RDATA_SEXPTYPE_PROMISE             5
#define RDATA_SEXPTYPE_LANGUAGE_OBJECT     6
#define RDATA_SEXPTYPE_SPECIAL_FUNCTION    7
#define RDATA_SEXPTYPE_BUILTIN_FUNCTION    8
#define RDATA_SEXPTYPE_CHARACTER_STRING    9
#define RDATA_SEXPTYPE_LOGICAL_VECTOR     10
#define RDATA_SEXPTYPE_INTEGER_VECTOR     13
#define RDATA_SEXPTYPE_REAL_VECTOR        14
#define RDATA_SEXPTYPE_COMPLEX_VECTOR     15
#define RDATA_SEXPTYPE_CHARACTER_VECTOR   16
#define RDATA_SEXPTYPE_DOT_DOT_DOT        17
#define RDATA_SEXPTYPE_ANY                18
#define RDATA_SEXPTYPE_GENERIC_VECTOR     19
#define RDATA_SEXPTYPE_EXPRESSION_VECTOR  20
#define RDATA_SEXPTYPE_BYTE_CODE          21
#define RDATA_SEXPTYPE_EXTERNAL_POINTER   22
#define RDATA_SEXPTYPE_WEAK_REFERENCE     23
#define RDATA_SEXPTYPE_RAW_VECTOR         24
#define RDATA_SEXPTYPE_S4_CLASS           25

#define RDATA_SEXPTYPE_FUN                99

#define RDATA_PSEUDO_SXP_REF                   255
#define RDATA_PSEUDO_SXP_NIL                   254
#define RDATA_PSEUDO_SXP_GLOBAL_ENVIRONMENT    253
#define RDATA_PSEUDO_SXP_UNBOUND_VALUE         252
#define RDATA_PSEUDO_SXP_MISSING_ARGUMENT      251
#define RDATA_PSEUDO_SXP_BASE_NAMESPACE        250
#define RDATA_PSEUDO_SXP_NAMESPACE             249
#define RDATA_PSEUDO_SXP_PACKAGE               248
#define RDATA_PSEUDO_SXP_PERSIST               247
#define RDATA_PSEUDO_SXP_CLASS_REF             246
#define RDATA_PSEUDO_SXP_GENERIC_REF           245
#define RDATA_PSEUDO_SXP_BYTE_CODE_REP_DEF     244
#define RDATA_PSEUDO_SXP_BYTE_CODE_REP_REF     243
#define RDATA_PSEUDO_SXP_EMPTY_ENVIRONMENT     242
#define RDATA_PSEUDO_SXP_BASE_ENVIRONMENT      241

#define RDATA_SEXPTYPE_LANGUAGE_OBJECT_ATTR    240
#define RDATA_SEXPTYPE_PAIRLIST_ATTR           239
#define RDATA_PSEUDO_SXP_ALTREP                238

#endif  // PANDAS_IO_RDATA_LIBRDATA_RDATA_INTERNAL_H_
