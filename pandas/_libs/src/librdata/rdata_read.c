/*
Copyright (c) 2020 Evan Miller
*/

//
//  rdata_rdata.c
//

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#if defined(_WIN32)
#include "win_iconv.h"
#elif __linux__
#include "unix_iconv.h"
#else
#include "<iconv.h>"
#endif

#include <errno.h>
#include <stdbool.h>

#if HAVE_BZIP2
#include <bzlib.h>
#endif

#if HAVE_APPLE_COMPRESSION
#include <compression.h>
#endif

#if HAVE_ZLIB
#include <zlib.h>
#endif

#if HAVE_LZMA
#include <lzma.h>
#endif

#include "rdata.h"
#include "rdata_internal.h"

#define RDATA_CLASS_POSIXCT 0x01
#define RDATA_CLASS_DATE    0x02

#define STREAM_BUFFER_SIZE   65536
#define MAX_ARRAY_DIMENSIONS     3

/* ICONV_CONST defined by autotools during configure according
 * to the current platform. Some people copy-paste the source code, so
 * provide some fallback logic */
#ifndef ICONV_CONST
#define ICONV_CONST
#endif

typedef struct rdata_atom_table_s {
    int    count;
    char **data;
} rdata_atom_table_t;

typedef struct rdata_ctx_s {
    int                          machine_needs_byteswap;
    rdata_table_handler          table_handler;
    rdata_column_handler         column_handler;
    rdata_column_name_handler    column_name_handler;
    rdata_column_name_handler    row_name_handler;
    rdata_text_value_handler     text_value_handler;
    rdata_text_value_handler     value_label_handler;
    rdata_column_handler         dim_handler;
    rdata_text_value_handler     dim_name_handler;
    rdata_error_handler       error_handler;
    void                        *user_ctx;
#if HAVE_BZIP2
    bz_stream                   *bz_strm;
#endif
#if HAVE_APPLE_COMPRESSION
    compression_stream          *compression_strm;
#endif
#if HAVE_ZLIB
    z_stream                    *z_strm;
#endif
#if HAVE_LZMA
    lzma_stream                 *lzma_strm;
#endif
    void                        *strm_buffer;
    rdata_io_t               *io;
    size_t                       bytes_read;

    rdata_atom_table_t          *atom_table;
    unsigned int                 column_class;

    iconv_t                      converter;

    int32_t                      dims[MAX_ARRAY_DIMENSIONS];
    bool                         is_dimnames;
} rdata_ctx_t;

static int atom_table_add(rdata_atom_table_t *table, char *key);
static char *atom_table_lookup(rdata_atom_table_t *table, int index);

static rdata_error_t read_environment(
    const char *table_name,
    rdata_ctx_t *ctx);
static rdata_error_t read_toplevel_object(
    const char *table_name,
    const char *key,
    rdata_ctx_t *ctx);
static rdata_error_t read_sexptype_header(
    rdata_sexptype_info_t *header,
    rdata_ctx_t *ctx);
static rdata_error_t read_length(
    int32_t *outLength,
    rdata_ctx_t *ctx);
static rdata_error_t read_string_vector_n(
    int attributes,
    int32_t length,
    rdata_text_value_handler text_value_handler,
    void *callback_ctx,
    rdata_ctx_t *ctx);
static rdata_error_t read_string_vector(
    int attributes,
    rdata_text_value_handler text_value_handler,
    void *callback_ctx,
    rdata_ctx_t *ctx);
static rdata_error_t read_value_vector(
    rdata_sexptype_header_t header,
    const char *name,
    rdata_ctx_t *ctx);
static rdata_error_t read_value_vector_cb(
    rdata_sexptype_header_t header,
    const char *name,
    rdata_column_handler column_handler,
    void *user_ctx,
    rdata_ctx_t *ctx);
static rdata_error_t read_character_string(
    char **key,
    rdata_ctx_t *ctx);
static rdata_error_t read_generic_list(
    int attributes,
    rdata_ctx_t *ctx);
static rdata_error_t read_altrep_vector(
    const char *name,
    rdata_ctx_t *ctx);
static rdata_error_t read_attributes(int (*handle_attribute)(
    char *key,
    rdata_sexptype_info_t val_info,
    rdata_ctx_t *ctx),
    rdata_ctx_t *ctx);
static rdata_error_t recursive_discard(
    rdata_sexptype_header_t sexptype_header,
    rdata_ctx_t *ctx);

static void *rdata_malloc(size_t len) {
    if (len == 0)
        return NULL;

    return malloc(len);
}

static void *rdata_realloc(void *buf, size_t len) {
    if (len == 0)
        return NULL;

    return realloc(buf, len);
}

static int atom_table_add(rdata_atom_table_t *table, char *key) {
    table->data = realloc(table->data, sizeof(char *) * (table->count + 1));
    table->data[table->count++] = strdup(key);
    return table->count;
}

static char *atom_table_lookup(rdata_atom_table_t *table, int index) {
    if (index <= 0 || index > table->count) {
        return NULL;
    }
    return table->data[(index-1)];
}

#if HAVE_BZIP2
static ssize_t read_st_bzip2(rdata_ctx_t *ctx, void *buffer, size_t len) {
    ssize_t bytes_written = 0;
    int error = 0;
    int result = BZ_OK;
    while (1) {
        ssize_t start_out = ctx->bz_strm->total_out_lo32 +
            ((ssize_t)ctx->bz_strm->total_out_hi32 << 32LL);

        ctx->bz_strm->next_out = (char *)buffer + bytes_written;
        ctx->bz_strm->avail_out = len - bytes_written;

        result = BZ2_bzDecompress(ctx->bz_strm);

        if (result != BZ_OK && result != BZ_STREAM_END) {
            error = -1;
            break;
        }

        bytes_written += ctx->bz_strm->total_out_lo32 +
            ((ssize_t)ctx->bz_strm->total_out_hi32 << 32LL) - start_out;

        if (result == BZ_STREAM_END)
            break;

        if (ctx->bz_strm->avail_in == 0) {
            int bytes_read = 0;
            bytes_read = ctx->io->read(
                ctx->strm_buffer,
                STREAM_BUFFER_SIZE,
                ctx->io->io_ctx);
            if (bytes_read < 0) {
                error = bytes_read;
                break;
            }
            if (bytes_read == 0)
                break;

            ctx->bz_strm->next_in = ctx->strm_buffer;
            ctx->bz_strm->avail_in = bytes_read;
        }
        if (bytes_written == len)
            break;
    }

    if (error != 0)
        return error;

    return bytes_written;
}
#endif /* HAVE_BZIP2 */

#if HAVE_APPLE_COMPRESSION
static ssize_t read_st_compression(
    rdata_ctx_t *ctx,
    void *buffer,
    size_t len
) {
    ssize_t bytes_written = 0;
    int error = 0;
    compression_status result = COMPRESSION_STATUS_OK;
    size_t start_size = len;

    ctx->compression_strm->dst_ptr = (unsigned char *)buffer;
    ctx->compression_strm->dst_size = len;

    while (1) {
        start_size = ctx->compression_strm->dst_size;

        result = compression_stream_process(ctx->compression_strm, 0);

        if (result == COMPRESSION_STATUS_OK) {
            bytes_written += start_size - ctx->compression_strm->dst_size;
        } else {
            error = -1;
            break;
        }

        if (ctx->compression_strm->src_size == 0) {
            int bytes_read = 0;
            bytes_read = ctx->io->read(
                ctx->compression_strm,
                STREAM_BUFFER_SIZE,
                ctx->io->io_ctx);
            if (bytes_read < 0) {
                error = bytes_read;
                break;
            }
            if (bytes_read == 0) {
                start_size = ctx->compression_strm->dst_size;
                result = compression_stream_process(
                    ctx->compression_strm,
                    COMPRESSION_STREAM_FINALIZE);
                if (result == COMPRESSION_STATUS_END) {
                    bytes_written += (
                        start_size - ctx->compression_strm->dst_size);
                } else {
                    error = -1;
                }
                break;
            }

            ctx->compression_strm->src_ptr = ctx->strm_buffer;
            ctx->compression_strm->src_size = bytes_read;
        }
        if (bytes_written == len)
            break;
    }

    if (error != 0)
        return error;

    return bytes_written;
}
#endif /* HAVE_APPLE_COMPRESSION */

#if HAVE_ZLIB
static ssize_t read_st_z(rdata_ctx_t *ctx, void *buffer, size_t len) {
    ssize_t bytes_written = 0;
    int error = 0;
    int result = Z_OK;
    while (1) {
        long start_out = ctx->z_strm->total_out;

        ctx->z_strm->next_out = (unsigned char *)buffer + bytes_written;
        ctx->z_strm->avail_out = len - bytes_written;

        result = inflate(ctx->z_strm, Z_SYNC_FLUSH);

        if (result != Z_OK && result != Z_STREAM_END) {
            error = -1;
            break;
        }

        bytes_written += ctx->z_strm->total_out - start_out;

        if (result == Z_STREAM_END)
            break;

        if (ctx->z_strm->avail_in == 0) {
            int bytes_read = 0;
            bytes_read = ctx->io->read(
                ctx->strm_buffer,
                STREAM_BUFFER_SIZE,
                ctx->io->io_ctx);
            if (bytes_read < 0) {
                error = bytes_read;
                break;
            }
            if (bytes_read == 0)
                break;

            ctx->z_strm->next_in = ctx->strm_buffer;
            ctx->z_strm->avail_in = bytes_read;
        }
        if (bytes_written == len)
            break;
    }

    if (error != 0)
        return error;

    return bytes_written;
}
#endif /* HAVE_ZLIB */

#if HAVE_LZMA
static ssize_t read_st_lzma(rdata_ctx_t *ctx, void *buffer, size_t len) {
    ssize_t bytes_written = 0;
    int error = 0;
    int result = LZMA_OK;
    while (1) {
        long start_out = ctx->lzma_strm->total_out;

        ctx->lzma_strm->next_out = (unsigned char *)buffer + bytes_written;
        ctx->lzma_strm->avail_out = len - bytes_written;

        result = lzma_code(ctx->lzma_strm, LZMA_RUN);

        if (result != LZMA_OK && result != LZMA_STREAM_END) {
            error = -1;
            break;
        }

        bytes_written += ctx->lzma_strm->total_out - start_out;

        if (result == LZMA_STREAM_END)
            break;

        if (ctx->lzma_strm->avail_in == 0) {
            int bytes_read = 0;
            bytes_read = ctx->io->read(
                ctx->strm_buffer,
                STREAM_BUFFER_SIZE,
                ctx->io->io_ctx);
            if (bytes_read < 0) {
                error = bytes_read;
                break;
            }
            if (bytes_read == 0)
                break;

            ctx->lzma_strm->next_in = ctx->strm_buffer;
            ctx->lzma_strm->avail_in = bytes_read;
        }
        if (bytes_written == len)
            break;
    }

    if (error != 0)
        return error;

    return bytes_written;
}
#endif /* HAVE_LZMA */

static ssize_t read_st(rdata_ctx_t *ctx, void *buffer, size_t len) {
    ssize_t bytes_read = 0;

    if (len == 0)
        return 0;

#if HAVE_BZIP2
    if (ctx->bz_strm) {
        bytes_read = read_st_bzip2(ctx, buffer, len);
    } else  // NOLINT
#endif
#if HAVE_APPLE_COMPRESSION
    if (ctx->compression_strm) {
        bytes_read = read_st_compression(ctx, buffer, len);
    } else  // NOLINT
#endif
#if HAVE_ZLIB
    if (ctx->z_strm) {
        bytes_read = read_st_z(ctx, buffer, len);
    } else  // NOLINT
#endif
#if HAVE_LZMA
    if (ctx->lzma_strm) {
        bytes_read = read_st_lzma(ctx, buffer, len);
    } else  // NOLINT
#endif
    {
        bytes_read = ctx->io->read(buffer, len, ctx->io->io_ctx);
    }

    if (bytes_read > 0) {
        ctx->bytes_read += bytes_read;
    }

    return bytes_read;
}

static int lseek_st(rdata_ctx_t *ctx, size_t len) {
    if (0
#if HAVE_BZIP2
            || ctx->bz_strm
#endif
#if HAVE_APPLE_COMPRESSION
            || ctx->compression_strm
#endif
#if HAVE_ZLIB
            || ctx->z_strm
#endif
#if HAVE_LZMA
            || ctx->lzma_strm
#endif
            ) {
        int retval = 0;
        char *buf = rdata_malloc(len);

        int result_st = read_st(ctx, buf, len);

        if (result_st > 0) {
            if (buf == NULL) {
                retval = -1;
            } else if ((size_t)result_st != len) {
                retval = -1;
            }
        } else {
            if (buf == NULL) {
                retval = -1;
            } else {
                retval = -1;
            }
        }

        if (buf)
            free(buf);

        return retval;
    }

    return ctx->io->seek(len, SEEK_CUR, ctx->io->io_ctx);
}

static rdata_error_t init_bz_stream(rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    ctx->strm_buffer = malloc(STREAM_BUFFER_SIZE);
    int bytes_read = ctx->io->read(
        ctx->strm_buffer,
        STREAM_BUFFER_SIZE,
        ctx->io->io_ctx);
    if (bytes_read <= 0) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

#if HAVE_BZIP2
    ctx->bz_strm = calloc(1, sizeof(bz_stream));
    ctx->bz_strm->next_in = ctx->strm_buffer;
    ctx->bz_strm->avail_in = bytes_read;

    if (BZ2_bzDecompressInit(ctx->bz_strm, 0, 0) != BZ_OK) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }
#else
    retval = RDATA_ERROR_UNSUPPORTED_COMPRESSION;
    goto cleanup;
#endif

cleanup:
    return retval;
}

static rdata_error_t init_z_stream(rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    ctx->strm_buffer = malloc(STREAM_BUFFER_SIZE);
    int bytes_read = ctx->io->read(
        ctx->strm_buffer,
        STREAM_BUFFER_SIZE,
        ctx->io->io_ctx);
    if (bytes_read <= 0) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

#if HAVE_ZLIB
    ctx->z_strm = calloc(1, sizeof(z_stream));
    ctx->z_strm->next_in = ctx->strm_buffer;
    ctx->z_strm->avail_in = bytes_read;

    if (inflateInit2(ctx->z_strm, (15+32)) != Z_OK) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }
#else
    retval = RDATA_ERROR_UNSUPPORTED_COMPRESSION;
    goto cleanup;
#endif

cleanup:
    return retval;
}

static rdata_error_t init_lzma_stream(rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    ctx->strm_buffer = malloc(STREAM_BUFFER_SIZE);
    int bytes_read = ctx->io->read(
        ctx->strm_buffer,
        STREAM_BUFFER_SIZE,
        ctx->io->io_ctx);
    if (bytes_read <= 0) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

#if HAVE_APPLE_COMPRESSION
    ctx->compression_strm = calloc(1, sizeof(compression_stream));

    if (compression_stream_init(
        ctx->compression_strm,
        COMPRESSION_STREAM_DECODE,
        COMPRESSION_LZMA) == COMPRESSION_STATUS_ERROR
    ) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }

    ctx->compression_strm->src_ptr = ctx->strm_buffer;
    ctx->compression_strm->src_size = bytes_read;
#elif HAVE_LZMA
    ctx->lzma_strm = calloc(1, sizeof(lzma_stream));

    if (lzma_stream_decoder(ctx->lzma_strm, UINT64_MAX, 0) != LZMA_OK) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }

    ctx->lzma_strm->next_in = ctx->strm_buffer;
    ctx->lzma_strm->avail_in = bytes_read;
#else
    retval = RDATA_ERROR_UNSUPPORTED_COMPRESSION;
    goto cleanup;
#endif

cleanup:
    return retval;
}

static rdata_error_t init_stream(rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    char header[5];

    if (ctx->io->read(
        &header,
        sizeof(header),
        ctx->io->io_ctx) != sizeof(header)
    ) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

    if (ctx->io->seek(0, SEEK_SET, ctx->io->io_ctx) == -1) {
        retval = RDATA_ERROR_SEEK;
        goto cleanup;
    }

    if (header[0] == 'B' && header[1] == 'Z' && header[2] == 'h' &&
            header[3] >= '0' && header[3] <= '9') {
        return init_bz_stream(ctx);
    }
    if (header[0] == '\x1f' && header[1] == '\x8b') {
        return init_z_stream(ctx);
    }
    if (strncmp("\xFD" "7zXZ", header, sizeof(header)) == 0) {
        return init_lzma_stream(ctx);
    }

cleanup:
    return retval;
}

static rdata_error_t reset_stream(rdata_ctx_t *ctx) {
#if HAVE_BZIP2
    if (ctx->bz_strm) {
        BZ2_bzDecompressEnd(ctx->bz_strm);
        free(ctx->bz_strm);
        ctx->bz_strm = NULL;
    }
#endif
#if HAVE_APPLE_COMPRESSION
    if (ctx->compression_strm) {
        compression_stream_destroy(ctx->compression_strm);
        free(ctx->compression_strm);
        ctx->compression_strm = NULL;
    }
#endif
#if HAVE_ZLIB
    if (ctx->z_strm) {
        inflateEnd(ctx->z_strm);
        free(ctx->z_strm);
        ctx->z_strm = NULL;
    }
#endif
#if HAVE_LZMA
    if (ctx->lzma_strm) {
        lzma_end(ctx->lzma_strm);
        free(ctx->lzma_strm);
        ctx->lzma_strm = NULL;
    }
#endif

    if (ctx->io->seek(0, SEEK_SET, ctx->io->io_ctx) == -1) {
        return RDATA_ERROR_SEEK;
    }
    return init_stream(ctx);
}

static rdata_error_t rdata_convert(
    char *dst,
    size_t dst_len,
    const char *src,
    size_t src_len,
    iconv_t converter
) {
    if (dst_len == 0) {
        return RDATA_ERROR_CONVERT_LONG_STRING;
    } else if (converter) {
        size_t dst_left = dst_len - 1;
        char *dst_end = dst;
        size_t status = iconv(converter, (
            ICONV_CONST char **)&src,
            &src_len,
            &dst_end,
            &dst_left);
        if (status == (size_t)-1) {
            if (errno == E2BIG) {
                return RDATA_ERROR_CONVERT_LONG_STRING;
            } else if (errno == EILSEQ) {
                return RDATA_ERROR_CONVERT_BAD_STRING;
            } else if (errno != EINVAL) {
                /* EINVAL indicates improper truncation; accept it */
                return RDATA_ERROR_CONVERT;
            }
        }
        dst[dst_len - dst_left - 1] = '\0';
    } else if (src_len + 1 > dst_len) {
        return RDATA_ERROR_CONVERT_LONG_STRING;
    } else {
        memcpy(dst, src, src_len);
        dst[src_len] = '\0';
    }
    return RDATA_OK;
}

rdata_ctx_t *rdata_ctx_init(rdata_io_t *io, const char *filename) {
    int fd = io->open(filename, io->io_ctx);
    if (fd == -1) {
        return NULL;
    }
    rdata_ctx_t *ctx = calloc(1, sizeof(rdata_ctx_t));
    rdata_atom_table_t *atom_table = malloc(sizeof(rdata_atom_table_t));

    atom_table->count = 0;
    atom_table->data = NULL;

    ctx->atom_table = atom_table;

    ctx->machine_needs_byteswap = 0;
    if (machine_is_little_endian()) {
        ctx->machine_needs_byteswap = 1;
    }

    ctx->io = io;

    return ctx;
}

void free_rdata_ctx(rdata_ctx_t *ctx) {
    if (ctx->io) {
        ctx->io->close(ctx->io->io_ctx);
    }
    if (ctx->atom_table) {
        if (ctx->atom_table->data) {
            int i;
            for (i=0; i < ctx->atom_table->count; i++)
                free(ctx->atom_table->data[i]);
            free(ctx->atom_table->data);
        }
        free(ctx->atom_table);
    }
#if HAVE_BZIP2
    if (ctx->bz_strm) {
        BZ2_bzDecompressEnd(ctx->bz_strm);
        free(ctx->bz_strm);
    }
#endif
#if HAVE_APPLE_COMPRESSION
    if (ctx->compression_strm) {
        compression_stream_destroy(ctx->compression_strm);
        free(ctx->compression_strm);
    }
#endif
#if HAVE_ZLIB
    if (ctx->z_strm) {
        inflateEnd(ctx->z_strm);
        free(ctx->z_strm);
    }
#endif
#if HAVE_LZMA
    if (ctx->lzma_strm) {
        lzma_end(ctx->lzma_strm);
        free(ctx->lzma_strm);
    }
#endif
    if (ctx->strm_buffer) {
        free(ctx->strm_buffer);
    }
    if (ctx->converter) {
        iconv_close(ctx->converter);
    }
    free(ctx);
}

rdata_error_t rdata_parse(
    rdata_parser_t *parser,
    const char *filename,
    void *user_ctx
) {
    int is_rdata = 0;
    rdata_error_t retval = RDATA_OK;
    rdata_v2_header_t v2_header;
    rdata_ctx_t *ctx = rdata_ctx_init(parser->io, filename);
    char *encoding = NULL;

    if (ctx == NULL) {
        retval = RDATA_ERROR_OPEN;
        goto cleanup;
    }

    ctx->user_ctx = user_ctx;
    ctx->table_handler = parser->table_handler;
    ctx->column_handler = parser->column_handler;
    ctx->column_name_handler = parser->column_name_handler;
    ctx->row_name_handler = parser->row_name_handler;
    ctx->text_value_handler = parser->text_value_handler;
    ctx->value_label_handler = parser->value_label_handler;
    ctx->dim_handler = parser->dim_handler;
    ctx->dim_name_handler = parser->dim_name_handler;
    ctx->error_handler = parser->error_handler;

    ctx->is_dimnames = false;

    if ((retval = init_stream(ctx)) != RDATA_OK) {
        goto cleanup;
    }

    char header_line[5];
    if (read_st(
        ctx, &header_line,
        sizeof(header_line)) != sizeof(header_line)
    ) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }
    if (memcmp("RDX", header_line, 3) == 0 && header_line[4] == '\n') {
        is_rdata = 1;
    } else {
        reset_stream(ctx);
    }

    if (read_st(ctx, &v2_header, sizeof(v2_header)) != sizeof(v2_header)) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

    if (ctx->machine_needs_byteswap) {
        v2_header.format_version = byteswap4(v2_header.format_version);
        v2_header.writer_version = byteswap4(v2_header.writer_version);
        v2_header.reader_version = byteswap4(v2_header.reader_version);
    }

    int32_t hdr_result = header_line[3] - '0';

    if (hdr_result > 0) {
        if (is_rdata && v2_header.format_version != (uint32_t)hdr_result) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }
    } else {
        if (is_rdata) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }
    }

    if (v2_header.format_version == 3) {
        retval = read_character_string(&encoding, ctx);
        if (retval != RDATA_OK)
            goto cleanup;

        if (strcmp("UTF-8", encoding) != 0) {
            if ((ctx->converter = iconv_open("UTF-8", encoding))
                == (iconv_t)-1
            ) {
                ctx->converter = NULL;
                retval = RDATA_ERROR_UNSUPPORTED_CHARSET;
                goto cleanup;
            }
        }
    }

    if (is_rdata) {
        retval = read_environment(NULL, ctx);
    } else {
        retval = read_toplevel_object(NULL, NULL, ctx);
    }
    if (retval != RDATA_OK)
        goto cleanup;

    char test;

    if (read_st(ctx, &test, 1) == 1) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }

cleanup:
    if (encoding)
        free(encoding);
    if (ctx) {
        free_rdata_ctx(ctx);
    }

    return retval;
}


static rdata_error_t read_toplevel_object(
    const char *table_name,
    const char *key,
    rdata_ctx_t *ctx
) {
    rdata_sexptype_info_t sexptype_info;
    rdata_error_t retval = RDATA_OK;

    sexptype_info.attributes = 0;
    sexptype_info.tag = 0;
    sexptype_info.ref = 0;
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;

    if (sexptype_info.header.type == RDATA_SEXPTYPE_REAL_VECTOR ||
            sexptype_info.header.type == RDATA_SEXPTYPE_INTEGER_VECTOR ||
            sexptype_info.header.type == RDATA_SEXPTYPE_LOGICAL_VECTOR) {
        if (table_name == NULL && ctx->table_handler) {
            if (ctx->table_handler(key, ctx->user_ctx)) {
                retval = RDATA_ERROR_USER_ABORT;
                goto cleanup;
            }
        }

        if ((retval = read_value_vector(
            sexptype_info.header,
            key,
            ctx)) != RDATA_OK
        )
            goto cleanup;
    } else if (sexptype_info.header.type == RDATA_SEXPTYPE_CHARACTER_VECTOR) {
        if (table_name == NULL && ctx->table_handler) {
            if (ctx->table_handler(key, ctx->user_ctx)) {
                retval = RDATA_ERROR_USER_ABORT;
                goto cleanup;
            }
        }
        int32_t length;

        if ((retval = read_length(&length, ctx)) != RDATA_OK)
            goto cleanup;

        if (ctx->column_handler) {
            if (ctx->column_handler(
                key,
                RDATA_TYPE_STRING, NULL,
                length, ctx->user_ctx)
            ) {
                retval = RDATA_ERROR_USER_ABORT;
                goto cleanup;
            }
        }

        if ((retval = read_string_vector_n(
            sexptype_info.header.attributes,
            length,
            ctx->text_value_handler,
            ctx->user_ctx, ctx)) != RDATA_OK)
            goto cleanup;
    } else if (sexptype_info.header.type == RDATA_PSEUDO_SXP_ALTREP) {
        if (table_name == NULL && ctx->table_handler) {
            if (ctx->table_handler(key, ctx->user_ctx)) {
                retval = RDATA_ERROR_USER_ABORT;
                goto cleanup;
            }
        }
        if ((retval = read_altrep_vector(key, ctx)) != RDATA_OK)
            goto cleanup;
    } else if (sexptype_info.header.type == RDATA_SEXPTYPE_GENERIC_VECTOR &&
            sexptype_info.header.object && sexptype_info.header.attributes) {
        if (table_name != NULL) {
            retval = recursive_discard(sexptype_info.header, ctx);
        } else {
            if (ctx->table_handler) {
                if (ctx->table_handler(key, ctx->user_ctx)) {
                    retval = RDATA_ERROR_USER_ABORT;
                    goto cleanup;
                }
            }
            retval = read_generic_list(sexptype_info.header.attributes, ctx);
        }
        if (retval != RDATA_OK)
            goto cleanup;
    } else {
        if ((retval = recursive_discard(sexptype_info.header, ctx))
            != RDATA_OK
        )
            goto cleanup;
    }

cleanup:

    return retval;
}

static rdata_error_t read_environment(
    const char *table_name,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    char *key = NULL;

    while (1) {
        rdata_sexptype_info_t sexptype_info;

        if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
            goto cleanup;

        if (sexptype_info.header.type == RDATA_PSEUDO_SXP_NIL)
            break;

        if (sexptype_info.header.type != RDATA_SEXPTYPE_PAIRLIST) {
            if ((retval = recursive_discard(
                sexptype_info.header,
                ctx)) != RDATA_OK)
                goto cleanup;
            continue;
        }

        if ((key = atom_table_lookup(
                ctx->atom_table,
                sexptype_info.ref)) == NULL) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }

        if ((retval = read_toplevel_object(table_name, key, ctx)) != RDATA_OK)
            goto cleanup;
    }

cleanup:

    return retval;
}

static rdata_error_t read_sexptype_header(
    rdata_sexptype_info_t *header_info,
    rdata_ctx_t *ctx
) {
    uint32_t sexptype;
    rdata_sexptype_header_t header;
    rdata_error_t retval = RDATA_OK;
    if (read_st(ctx, &sexptype, sizeof(sexptype)) != sizeof(sexptype)) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }
    if (ctx->machine_needs_byteswap)
        sexptype = byteswap4(sexptype);

    memcpy(&header, &sexptype, sizeof(sexptype));
    uint32_t attributes = 0, tag = 0, ref = 0;

    if (header.type == RDATA_SEXPTYPE_PAIRLIST_ATTR) {
        header.attributes = 1;
        header.type = RDATA_SEXPTYPE_PAIRLIST;
    }
    if (header.type == RDATA_SEXPTYPE_LANGUAGE_OBJECT_ATTR) {
        header.attributes = 1;
        header.type = RDATA_SEXPTYPE_LANGUAGE_OBJECT;
    }
    if (header.type == RDATA_SEXPTYPE_PAIRLIST) {
        if (header.attributes) {
            if (read_st(
                ctx,
                &attributes,
                sizeof(attributes)) != sizeof(attributes)
            ) {
                retval = RDATA_ERROR_READ;
                goto cleanup;
            }
            if (ctx->machine_needs_byteswap) {
                header_info->attributes = byteswap4(header_info->attributes);
            }
        }
        if (header.tag) {
            if (read_st(ctx, &tag, sizeof(tag)) != sizeof(tag)) {
                retval = RDATA_ERROR_READ;
                goto cleanup;
            }
            if (ctx->machine_needs_byteswap)
                tag = byteswap4(tag);
        }

        if (tag == 1) {
            rdata_sexptype_info_t key_info;

            if ((retval = read_sexptype_header(&key_info, ctx)) != RDATA_OK)
                goto cleanup;

            if (key_info.header.type != RDATA_SEXPTYPE_CHARACTER_STRING) {
                retval = RDATA_ERROR_PARSE;
                goto cleanup;
            }

            char *key = NULL;
            if ((retval = read_character_string(&key, ctx)) != RDATA_OK)
                goto cleanup;

            ref = atom_table_add(ctx->atom_table, key);

            free(key);
        } else if ((tag & 0xFF) == RDATA_PSEUDO_SXP_REF) {
            ref = (tag >> 8);
        }
    }
    if (header.type == RDATA_PSEUDO_SXP_REF) {
        ref = (sexptype >> 8);
    }

    header_info->header = header;
    header_info->attributes = attributes;
    header_info->tag = tag;
    header_info->ref = ref;

cleanup:

    return retval;
}

static int handle_class_name(const char *buf, int i, void *ctx) {
    unsigned int *column_class = (unsigned int *)ctx;
    if (buf) {
        if (strcmp(buf, "POSIXct") == 0) {
            *column_class |= RDATA_CLASS_POSIXCT;
        }
        if (strcmp(buf, "Date") == 0) {
            *column_class |= RDATA_CLASS_DATE;
        }
    }
    return RDATA_OK;
}

static int handle_vector_attribute(
    char *key,
    rdata_sexptype_info_t val_info,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    if (strcmp(key, "levels") == 0) {
        retval = read_string_vector(
            val_info.header.attributes,
            ctx->value_label_handler,
            ctx->user_ctx, ctx);
    } else if (strcmp(key, "class") == 0) {
        ctx->column_class = 0;
        retval = read_string_vector(
            val_info.header.attributes,
            &handle_class_name,
            &ctx->column_class, ctx);
    } else if (strcmp(key, "dim") == 0) {
        if (val_info.header.type == RDATA_SEXPTYPE_INTEGER_VECTOR) {
            int32_t length;
            if ((retval = read_length(&length, ctx)) != RDATA_OK)
                goto cleanup;

            if ((uint32_t)length <= sizeof(ctx->dims)/sizeof(ctx->dims[0])) {
                int buf_len = length * sizeof(int32_t);
                if (read_st(ctx, ctx->dims, buf_len) != buf_len) {
                    retval = RDATA_ERROR_READ;
                    goto cleanup;
                }
                if (ctx->machine_needs_byteswap) {
                    int i;
                    for (i=0; i < length; i++) {
                        ctx->dims[i] = byteswap4(ctx->dims[i]);
                    }
                }
                if (ctx->dim_handler) {
                    if (ctx->dim_handler(
                        key,
                        RDATA_TYPE_INT32,
                        ctx->dims, length,
                        ctx->user_ctx)
                    ) {
                        retval = RDATA_ERROR_USER_ABORT;
                    }
                }
            }
        }
    } else if (strcmp(key, "dimnames") == 0) {
        ctx->is_dimnames = true;
        retval = read_generic_list(val_info.header.attributes, ctx);
    } else {
        retval = recursive_discard(val_info.header, ctx);
    }
cleanup:
    return retval;
}

static rdata_error_t read_character_string(char **key, rdata_ctx_t *ctx) {
    uint32_t length;
    char *string = NULL;
    char *utf8_string = NULL;
    rdata_error_t retval = RDATA_OK;

    if (read_st(ctx, &length, sizeof(length)) != sizeof(length)) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

    if (ctx->machine_needs_byteswap)
        length = byteswap4(length);

    if ((int32_t)length == -1 || length == 0) {
        *key = strdup("");
        return RDATA_OK;
    }

    if (length < 0) {
        return RDATA_ERROR_PARSE;
    }

    if ((string = rdata_malloc(length)) == NULL) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }

    if (read_st(ctx, string, length) != length) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

    if ((utf8_string = rdata_malloc(4*length+1)) == NULL) {
        retval = RDATA_ERROR_MALLOC;
        goto cleanup;
    }

    retval = rdata_convert(
        utf8_string,
        4 * length + 1,
        string, length,
        ctx->converter);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    if (string)
        free(string);

    if (retval == RDATA_OK) {
        *key = utf8_string;
    } else if (utf8_string) {
        free(utf8_string);
    }

    return retval;
}

static int handle_data_frame_attribute(
    char *key,
    rdata_sexptype_info_t val_info,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;

    if (strcmp(key, "names") == 0 &&
        val_info.header.type == RDATA_SEXPTYPE_CHARACTER_VECTOR
    ) {
        retval = read_string_vector(
            val_info.header.attributes,
            ctx->column_name_handler, ctx->user_ctx, ctx);
    } else if (strcmp(key, "row.names") == 0 &&
        val_info.header.type == RDATA_SEXPTYPE_CHARACTER_VECTOR
    ) {
        retval = read_string_vector(
        val_info.header.attributes,
        ctx->row_name_handler,
        ctx->user_ctx, ctx);
    } else if (strcmp(key, "label.table") == 0) {
        retval = recursive_discard(val_info.header, ctx);
    } else {
        retval = recursive_discard(val_info.header, ctx);
    }

    return retval;
}

static rdata_error_t read_attributes(int (*handle_attribute)(
    char *key,
    rdata_sexptype_info_t val_info,
    rdata_ctx_t *ctx),
    rdata_ctx_t *ctx
) {
    rdata_sexptype_info_t pairlist_info, val_info;
    rdata_error_t retval = RDATA_OK;
    char *key = NULL;

    retval = read_sexptype_header(&pairlist_info, ctx);
    if (retval != RDATA_OK)
        goto cleanup;

    while (pairlist_info.header.type == RDATA_SEXPTYPE_PAIRLIST) {
        /* value */
        if ((retval = read_sexptype_header(&val_info, ctx)) != RDATA_OK)
            goto cleanup;

        if (handle_attribute) {
            if ((key = atom_table_lookup(
                ctx->atom_table, pairlist_info.ref)) == NULL) {
                retval = RDATA_ERROR_PARSE;
                goto cleanup;
            }
            if ((retval = handle_attribute(key, val_info, ctx)) != RDATA_OK)
                goto cleanup;
        } else {
            if ((retval = recursive_discard(
                    val_info.header,
                    ctx)) != RDATA_OK
                )
                goto cleanup;
        }

        /* next */
        if ((retval = read_sexptype_header(&pairlist_info, ctx)) != RDATA_OK)
            goto cleanup;
    }

cleanup:
    return retval;
}

static rdata_error_t read_wrap_real(const char *name, rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    rdata_sexptype_info_t sexptype_info;
    /* pairlist */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_SEXPTYPE_PAIRLIST) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }
    /* representation */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;

    if ((retval = read_value_vector(
            sexptype_info.header,
            name,
            ctx)) != RDATA_OK
        )
        goto cleanup;

    /* alt representation */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if ((retval = recursive_discard(sexptype_info.header, ctx)) != RDATA_OK)
        goto cleanup;

    /* nil */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_PSEUDO_SXP_NIL) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }

cleanup:
    return retval;
}

static rdata_error_t read_compact_intseq(
    const char *name,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    rdata_sexptype_info_t sexptype_info;
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;

    int32_t length;
    if ((retval = read_length(&length, ctx)) != RDATA_OK)
        goto cleanup;
    if (length != 3) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }

    double vals[3];
    if (read_st(ctx, vals, sizeof(vals)) != sizeof(vals)) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }
    if (ctx->machine_needs_byteswap) {
        vals[0] = byteswap_double(vals[0]);
        vals[1] = byteswap_double(vals[1]);
        vals[2] = byteswap_double(vals[2]);
    }

    if (sexptype_info.header.attributes) {
        if ((retval = read_attributes(
            &handle_vector_attribute, ctx)) != RDATA_OK
        )
            goto cleanup;
    }

    if (ctx->column_handler) {
        int32_t *integers = rdata_malloc(vals[0] * sizeof(int32_t));
        int32_t val = vals[1];
        for (int i=0; i < vals[0]; i++) {
            integers[i] = val;
            val += vals[2];
        }
        int cb_retval = ctx->column_handler(
            name,
            RDATA_TYPE_INT32,
            integers,
            vals[0], ctx->user_ctx);
        free(integers);
        if (cb_retval) {
            retval = RDATA_ERROR_USER_ABORT;
            goto cleanup;
        }
    }

    /* nil */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_PSEUDO_SXP_NIL) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }
cleanup:
    return retval;
}

static int deferred_string_handler(
    const char *name,
    enum rdata_type_e type,
    void *vals,
    long length,
    void *user_ctx
) {
    rdata_ctx_t *ctx = (rdata_ctx_t *)user_ctx;
    if (ctx->column_handler)
        ctx->column_handler(
            name,
            RDATA_TYPE_STRING,
            NULL,
            length,
            ctx->user_ctx);
    if (ctx->text_value_handler) {
        for (int i=0; i < length; i++) {
            char buf[128] = { 0 };
            if (type == RDATA_TYPE_REAL) {
                snprintf(buf, sizeof(buf), "%.0lf", ((double *)vals)[i]);
            } else if (type == RDATA_TYPE_INT32) {
                snprintf(buf, sizeof(buf), "%d", ((int32_t *)vals)[i]);
            }
            ctx->text_value_handler(buf, i, ctx->user_ctx);
        }
    }
    return 0;
}

static rdata_error_t read_deferred_string(
    const char *name,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    rdata_sexptype_info_t sexptype_info;
    /* pairlist */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_SEXPTYPE_PAIRLIST) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }
    /* representation */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;

    if ((retval = read_value_vector_cb(
        sexptype_info.header,
        name,
        &deferred_string_handler,
        ctx,
        ctx)) != RDATA_OK
    )
        goto cleanup;

    /* alt representation */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if ((retval = recursive_discard(sexptype_info.header, ctx)) != RDATA_OK)
        goto cleanup;

    /* nil */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_PSEUDO_SXP_NIL) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }

cleanup:
    return retval;
}

static rdata_error_t read_altrep_vector(
    const char *name,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    rdata_sexptype_info_t sexptype_info;
    /* pairlist */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_SEXPTYPE_PAIRLIST) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }
    /* class name */
    char *class = NULL;
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type == RDATA_SEXPTYPE_SYMBOL) {
        if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
            goto cleanup;
        if (sexptype_info.header.type != RDATA_SEXPTYPE_CHARACTER_STRING) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }
        if ((retval = read_character_string(&class, ctx)) != RDATA_OK)
            goto cleanup;

        atom_table_add(ctx->atom_table, class);
    } else if (sexptype_info.header.type == RDATA_PSEUDO_SXP_REF) {
        if ((class = atom_table_lookup(
            ctx->atom_table,
            sexptype_info.ref)) == NULL
        ) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }
    } else {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }

    /* package and class ID */
    if ((retval = read_sexptype_header(&sexptype_info, ctx)) != RDATA_OK)
        goto cleanup;
    if (sexptype_info.header.type != RDATA_SEXPTYPE_PAIRLIST) {
        retval = RDATA_ERROR_PARSE;
        goto cleanup;
    }
    if ((retval = recursive_discard(sexptype_info.header, ctx)) != RDATA_OK)
        goto cleanup;

    if (strcmp(class, "wrap_real") == 0) {
        if ((retval = read_wrap_real(name, ctx)) != RDATA_OK)
            goto cleanup;
    } else if (strcmp(class, "compact_intseq") == 0) {
        if ((retval = read_compact_intseq(name, ctx)) != RDATA_OK)
            goto cleanup;
    } else if (strcmp(class, "deferred_string") == 0) {
        if ((retval = read_deferred_string(name, ctx)) != RDATA_OK)
            goto cleanup;
    } else {
        if (ctx->error_handler) {
            char error_buf[1024];
            snprintf(
                error_buf,
                sizeof(error_buf),
                "Unrecognized ALTREP class: %s\n",
                class);
            ctx->error_handler(error_buf, ctx->user_ctx);
        }
        retval = RDATA_ERROR_UNSUPPORTED_STORAGE_CLASS;
    }
cleanup:
    return retval;
}

static rdata_error_t read_generic_list(int attributes, rdata_ctx_t *ctx) {
    rdata_error_t retval = RDATA_OK;
    int32_t length;
    unsigned int i;
    rdata_sexptype_info_t sexptype_info;

    if ((retval = read_length(&length, ctx)) != RDATA_OK)
        goto cleanup;

    for (i=0; i < (uint32_t)length; i++) {
        if ((retval = read_sexptype_header(
            &sexptype_info, ctx)) != RDATA_OK
        )
            goto cleanup;

        if (sexptype_info.header.type == RDATA_SEXPTYPE_CHARACTER_VECTOR) {
            int32_t vec_length;

            if ((retval = read_length(&vec_length, ctx)) != RDATA_OK)
                goto cleanup;
            if (ctx->is_dimnames) {
                retval = read_string_vector_n(
                    sexptype_info.header.attributes,
                    vec_length,
                    ctx->dim_name_handler,
                    ctx->user_ctx, ctx);
            } else {
                if (ctx->column_handler) {
                    if (ctx->column_handler(
                        NULL,
                        RDATA_TYPE_STRING,
                        NULL,
                        vec_length,
                        ctx->user_ctx)
                    ) {
                        retval = RDATA_ERROR_USER_ABORT;
                        goto cleanup;
                    }
                }
                retval = read_string_vector_n(
                    sexptype_info.header.attributes,
                    vec_length,
                    ctx->text_value_handler,
                    ctx->user_ctx, ctx);
            }
        } else if (sexptype_info.header.type == RDATA_PSEUDO_SXP_ALTREP) {
            retval = read_altrep_vector(NULL, ctx);
        } else if (sexptype_info.header.type == RDATA_PSEUDO_SXP_NIL) {
            if (ctx->is_dimnames &&
                ctx->dim_name_handler &&
                i < sizeof(ctx->dims)/sizeof(ctx->dims[0])
            ) {
                int j;
                for (j=0; j < ctx->dims[i]; j++) {
                    ctx->dim_name_handler(NULL, j, ctx->user_ctx);
                }
            }
        } else {
            retval = read_value_vector(sexptype_info.header, NULL, ctx);
        }
        if (retval != RDATA_OK)
            goto cleanup;
    }

    if (attributes) {
        if ((retval = read_attributes(
            &handle_data_frame_attribute,
            ctx)) != RDATA_OK
        )
            goto cleanup;
    }

cleanup:

    if (ctx->is_dimnames)
        ctx->is_dimnames = false;

    return retval;
}

static rdata_error_t read_length(int32_t *outLength, rdata_ctx_t *ctx) {
    int32_t length;
    rdata_error_t retval = RDATA_OK;

    if (read_st(ctx, &length, sizeof(length)) != sizeof(length)) {
        retval = RDATA_ERROR_READ;
        goto cleanup;
    }

    if (ctx->machine_needs_byteswap)
        length = byteswap4(length);

    if (outLength)
        *outLength = length;

cleanup:

    return retval;
}

static rdata_error_t read_string_vector_n(
    int attributes,
    int32_t length,
    rdata_text_value_handler text_value_handler,
    void *callback_ctx,
    rdata_ctx_t *ctx
) {
    int32_t string_length;
    rdata_error_t retval = RDATA_OK;
    rdata_sexptype_info_t info;
    size_t buffer_size = 4096;
    char *buffer = NULL;
    size_t utf8_buffer_size = 16384;
    char *utf8_buffer = NULL;
    int i;

    buffer = rdata_malloc(buffer_size);
    if (ctx->converter)
        utf8_buffer = rdata_malloc(utf8_buffer_size);

    for (i=0; i < length; i++) {
        if ((retval = read_sexptype_header(&info, ctx)) != RDATA_OK)
            goto cleanup;

        if (info.header.type != RDATA_SEXPTYPE_CHARACTER_STRING) {
            retval = RDATA_ERROR_PARSE;
            goto cleanup;
        }

        if ((retval = read_length(&string_length, ctx)) != RDATA_OK)
            goto cleanup;

        int32_t str_len_calc = string_length + 1;
        if (str_len_calc > 0) {
            if ((uint32_t)str_len_calc > buffer_size) {
                buffer_size = str_len_calc;
                if ((buffer = rdata_realloc(buffer, buffer_size)) == NULL) {
                    retval = RDATA_ERROR_MALLOC;
                    goto cleanup;
                }
            }
        }

        if (string_length >= 0) {
            if (read_st(ctx, buffer, string_length) != string_length) {
                retval = RDATA_ERROR_READ;
                goto cleanup;
            }
            buffer[string_length] = '\0';
        }

        if (text_value_handler) {
            int cb_retval = 0;
            if (string_length < 0) {
                cb_retval = text_value_handler(NULL, i, callback_ctx);
            } else if (!ctx->converter) {
                cb_retval = text_value_handler(buffer, i, callback_ctx);
            } else {
                int32_t str_len_calc = 4*string_length + 1;
                if (str_len_calc >= 0) {
                    if ((uint32_t)str_len_calc > utf8_buffer_size) {
                        utf8_buffer_size = str_len_calc;
                        if ((utf8_buffer = rdata_realloc(
                            utf8_buffer, utf8_buffer_size)) == NULL
                        ) {
                            retval = RDATA_ERROR_MALLOC;
                            goto cleanup;
                        }
                    }
                }

                retval = rdata_convert(
                    utf8_buffer,
                    utf8_buffer_size,
                    buffer, string_length,
                    ctx->converter);
                if (retval != RDATA_OK)
                    goto cleanup;

                cb_retval = text_value_handler(utf8_buffer, i, callback_ctx);
            }
            if (cb_retval) {
                retval = RDATA_ERROR_USER_ABORT;
                goto cleanup;
            }
        }
    }

    if (attributes) {
        if ((retval = read_attributes(
            &handle_vector_attribute,
            ctx)) != RDATA_OK)
            goto cleanup;
    }

cleanup:

    if (buffer)
        free(buffer);
    if (utf8_buffer)
        free(utf8_buffer);

    return retval;
}

static rdata_error_t read_string_vector(
    int attributes,
    rdata_text_value_handler text_value_handler,
    void *callback_ctx,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    int32_t length;

    if ((retval = read_length(&length, ctx)) != RDATA_OK)
        return retval;

    return read_string_vector_n(
        attributes,
        length,
        text_value_handler,
        callback_ctx,
        ctx);
}

static rdata_error_t read_value_vector_cb(
    rdata_sexptype_header_t header,
    const char *name,
    rdata_column_handler column_handler,
    void *user_ctx,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    int32_t length;
    size_t input_elem_size = 0;
    void *vals = NULL;
    size_t buf_len = 0;
    enum rdata_type_e output_data_type;
    unsigned int i;

    switch (header.type) {
        case RDATA_SEXPTYPE_REAL_VECTOR:
            input_elem_size = sizeof(double);
            output_data_type = RDATA_TYPE_REAL;
            break;
        case RDATA_SEXPTYPE_INTEGER_VECTOR:
            input_elem_size = sizeof(int32_t);
            output_data_type = RDATA_TYPE_INT32;
            break;
        case RDATA_SEXPTYPE_LOGICAL_VECTOR:
            input_elem_size = sizeof(int32_t);
            output_data_type = RDATA_TYPE_LOGICAL;
            break;
        default:
            retval = RDATA_ERROR_PARSE;
            break;
    }
    if (retval != RDATA_OK)
        goto cleanup;

    if ((retval = read_length(&length, ctx)) != RDATA_OK)
        goto cleanup;

    buf_len = length * input_elem_size;

    if (buf_len) {
        vals = rdata_malloc(buf_len);
        if (vals == NULL) {
            retval = RDATA_ERROR_MALLOC;
            goto cleanup;
        }

        ssize_t result_st = read_st(ctx, vals, buf_len);

        if (result_st > 0) {
            if ((size_t)result_st != buf_len) {
                retval = RDATA_ERROR_READ;
                goto cleanup;
            }
        } else {
            retval = RDATA_ERROR_READ;
            goto cleanup;
        }

        if (ctx->machine_needs_byteswap) {
            if (input_elem_size == sizeof(double)) {
                double *d_vals = (double *)vals;
                for (i=0; i < buf_len/sizeof(double); i++) {
                    d_vals[i] = byteswap_double(d_vals[i]);
                }
            } else {
                uint32_t *i_vals = (uint32_t *)vals;
                for (i=0; i < buf_len/sizeof(uint32_t); i++) {
                    i_vals[i] = byteswap4(i_vals[i]);
                }
            }
        }
    }

    ctx->column_class = 0;
    if (header.attributes) {
        if ((retval = read_attributes(
               &handle_vector_attribute,
               ctx)) != RDATA_OK)
            goto cleanup;
    }
    if (ctx->column_class == RDATA_CLASS_POSIXCT)
        output_data_type = RDATA_TYPE_TIMESTAMP;
    if (ctx->column_class == RDATA_CLASS_DATE)
        output_data_type = RDATA_TYPE_DATE;

    if (column_handler) {
        if (column_handler(name, output_data_type, vals, length, user_ctx)) {
            retval = RDATA_ERROR_USER_ABORT;
            goto cleanup;
        }
    }

cleanup:
    if (vals)
        free(vals);

    return retval;
}

static rdata_error_t read_value_vector(
    rdata_sexptype_header_t header,
    const char *name,
    rdata_ctx_t *ctx
) {
    return read_value_vector_cb(
        header,
        name,
        ctx->column_handler,
        ctx->user_ctx, ctx);
}

static rdata_error_t discard_vector(
    rdata_sexptype_header_t sexptype_header,
    size_t element_size,
    rdata_ctx_t *ctx
) {
    int32_t length;
    rdata_error_t retval = RDATA_OK;

    if ((retval = read_length(&length, ctx)) != RDATA_OK)
        goto cleanup;

    if (length > 0) {
        if (lseek_st(ctx, length * element_size) == -1) {
            return RDATA_ERROR_SEEK;
        }
    } else if (ctx->error_handler) {
        char error_buf[1024];
        snprintf(
            error_buf,
            sizeof(error_buf),
            "Vector with non-positive length: %d\n",
            length);
        ctx->error_handler(error_buf, ctx->user_ctx);
    }

    if (sexptype_header.attributes) {
        rdata_sexptype_info_t temp_info;
        if ((retval = read_sexptype_header(&temp_info, ctx)) != RDATA_OK)
            goto cleanup;

        retval = recursive_discard(temp_info.header, ctx);
    }

cleanup:

    return retval;
}

static rdata_error_t discard_character_string(
    int add_to_table,
    rdata_ctx_t *ctx
) {
    rdata_error_t retval = RDATA_OK;
    char *key = NULL;

    if ((retval = read_character_string(&key, ctx)) != RDATA_OK)
        goto cleanup;

    if (strlen(key) > 0 && add_to_table) {
        atom_table_add(ctx->atom_table, key);
    }

    free(key);

cleanup:

    return retval;
}

static rdata_error_t discard_pairlist(
    rdata_sexptype_header_t sexptype_header,
    rdata_ctx_t *ctx
) {
    rdata_sexptype_info_t temp_info;
    rdata_error_t error = 0;
    while (1) {
        switch (sexptype_header.type) {
            case RDATA_SEXPTYPE_PAIRLIST:
                /* value */
                if ((error = read_sexptype_header(
                    &temp_info,
                    ctx)) != RDATA_OK)
                    return error;
                if ((error = recursive_discard(
                    temp_info.header,
                    ctx)) != RDATA_OK)
                    return error;

                /* tail */
                if ((error = read_sexptype_header(
                    &temp_info,
                    ctx)) != RDATA_OK)
                    return error;
                sexptype_header = temp_info.header;
                break;
            case RDATA_PSEUDO_SXP_NIL:
                goto done;
            default:
                return RDATA_ERROR_PARSE;
        }
    }
done:

    return 0;
}

static rdata_error_t recursive_discard(
    rdata_sexptype_header_t sexptype_header,
    rdata_ctx_t *ctx
) {
    uint32_t length;
    rdata_sexptype_info_t info;
    rdata_sexptype_info_t prot, tag;

    rdata_error_t error = 0;
    unsigned int i;

    switch (sexptype_header.type) {
        case RDATA_SEXPTYPE_SYMBOL:
            if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                goto cleanup;

            if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                goto cleanup;
            break;
        case RDATA_PSEUDO_SXP_PERSIST:
        case RDATA_PSEUDO_SXP_NAMESPACE:
        case RDATA_PSEUDO_SXP_PACKAGE:
            if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                goto cleanup;

            if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                goto cleanup;
            break;
        case RDATA_SEXPTYPE_BUILTIN_FUNCTION:
        case RDATA_SEXPTYPE_SPECIAL_FUNCTION:
            error = discard_character_string(0, ctx);
            break;
        case RDATA_SEXPTYPE_PAIRLIST:
            error = discard_pairlist(sexptype_header, ctx);
            break;
        case RDATA_SEXPTYPE_CHARACTER_STRING:
            error = discard_character_string(1, ctx);
            break;
        case RDATA_SEXPTYPE_RAW_VECTOR:
            error = discard_vector(sexptype_header, 1, ctx);
            break;
        case RDATA_SEXPTYPE_LOGICAL_VECTOR:
            error = discard_vector(sexptype_header, 4, ctx);
            break;
        case RDATA_SEXPTYPE_INTEGER_VECTOR:
            error = discard_vector(sexptype_header, 4, ctx);
            break;
        case RDATA_SEXPTYPE_REAL_VECTOR:
            error = discard_vector(sexptype_header, 8, ctx);
            break;
        case RDATA_SEXPTYPE_COMPLEX_VECTOR:
            error = discard_vector(sexptype_header, 16, ctx);
            break;
        case RDATA_SEXPTYPE_CHARACTER_VECTOR:
        case RDATA_SEXPTYPE_GENERIC_VECTOR:
        case RDATA_SEXPTYPE_EXPRESSION_VECTOR:
            if (read_st(ctx, &length, sizeof(length)) != sizeof(length)) {
                return RDATA_ERROR_READ;
            }
            if (ctx->machine_needs_byteswap)
                length = byteswap4(length);

            for (i=0; i < length; i++) {
                if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                    goto cleanup;

                if (sexptype_header.type == RDATA_SEXPTYPE_CHARACTER_VECTOR) {
                    if (info.header.type != RDATA_SEXPTYPE_CHARACTER_STRING) {
                        error = RDATA_ERROR_PARSE;
                        goto cleanup;
                    }

                    if ((error = discard_character_string(0, ctx)) != RDATA_OK)
                        goto cleanup;
                } else if ((error = recursive_discard(
                        info.header,
                        ctx)) != RDATA_OK) {
                    goto cleanup;
                }
            }
            if (sexptype_header.attributes) {
                if ((error = read_attributes(NULL, ctx)) != RDATA_OK)
                    goto cleanup;
            }
            break;
        case RDATA_SEXPTYPE_DOT_DOT_DOT:
        case RDATA_SEXPTYPE_PROMISE:
        case RDATA_SEXPTYPE_LANGUAGE_OBJECT:
        case RDATA_SEXPTYPE_CLOSURE:
            if (sexptype_header.attributes) {
                if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                    goto cleanup;

                if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                    goto cleanup;
            }
            if (sexptype_header.tag) {
                if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                    goto cleanup;

                if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                    goto cleanup;
            }
            /* CAR */
            if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                goto cleanup;

            if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                goto cleanup;

            /* CDR */
            if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                goto cleanup;

            if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                goto cleanup;
            break;
        case RDATA_SEXPTYPE_EXTERNAL_POINTER:
            read_sexptype_header(&prot, ctx);
            recursive_discard(prot.header, ctx);

            read_sexptype_header(&tag, ctx);
            recursive_discard(tag.header, ctx);
            break;
        case RDATA_SEXPTYPE_ENVIRONMENT:
            /* locked */
            if (lseek_st(ctx, sizeof(uint32_t)) == -1) {
                return RDATA_ERROR_SEEK;
            }

            rdata_sexptype_info_t enclosure, frame, hash_table, attributes;
            read_sexptype_header(&enclosure, ctx);
            recursive_discard(enclosure.header, ctx);

            read_sexptype_header(&frame, ctx);
            recursive_discard(frame.header, ctx);

            read_sexptype_header(&hash_table, ctx);
            recursive_discard(hash_table.header, ctx);

            read_sexptype_header(&attributes, ctx);
            recursive_discard(attributes.header, ctx);
            /*
             if (sexptype_header.attributes) {
             if (lseek(ctx->fd, sizeof(uint32_t), SEEK_CUR) == -1) {
             return RDATA_ERROR_SEEK;
             }
             } */
            break;
        case RDATA_PSEUDO_SXP_REF:
        case RDATA_PSEUDO_SXP_NIL:
        case RDATA_PSEUDO_SXP_GLOBAL_ENVIRONMENT:
        case RDATA_PSEUDO_SXP_UNBOUND_VALUE:
        case RDATA_PSEUDO_SXP_MISSING_ARGUMENT:
        case RDATA_PSEUDO_SXP_BASE_NAMESPACE:
        case RDATA_PSEUDO_SXP_EMPTY_ENVIRONMENT:
        case RDATA_PSEUDO_SXP_BASE_ENVIRONMENT:
            break;
        case RDATA_PSEUDO_SXP_ALTREP:
            /* class, package, type */
            if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                goto cleanup;
            if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                goto cleanup;

            while (1) {
                if ((error = read_sexptype_header(&info, ctx)) != RDATA_OK)
                    goto cleanup;
                if (info.header.type == RDATA_SEXPTYPE_PAIRLIST)
                    continue;
                if (info.header.type == RDATA_PSEUDO_SXP_NIL)
                    break;
                if ((error = recursive_discard(info.header, ctx)) != RDATA_OK)
                    goto cleanup;
            }
            break;
        default:
            if (ctx->error_handler) {
                char error_buf[1024];
                snprintf(
                    error_buf,
                    sizeof(error_buf),
                    "Unhandled S-Expression: %d",
                    sexptype_header.type);
                ctx->error_handler(error_buf, ctx->user_ctx);
            }
            return RDATA_ERROR_UNSUPPORTED_S_EXPRESSION;
    }
cleanup:

    return error;
}
