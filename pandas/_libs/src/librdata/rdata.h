/*
Copyright (c) 2020 Evan Miller
*/

#ifndef PANDAS__LIBS_SRC_LIBRDATA_RDATA_H_
#define PANDAS__LIBS_SRC_LIBRDATA_RDATA_H_

#include <stdint.h>
#include <sys/types.h>
#include <stdio.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum rdata_type_e {
    RDATA_TYPE_STRING,
    RDATA_TYPE_INT32,
    RDATA_TYPE_REAL,
    RDATA_TYPE_LOGICAL,
    RDATA_TYPE_TIMESTAMP,
    RDATA_TYPE_DATE
} rdata_type_t;

typedef enum rdata_error_e {
    RDATA_OK,
    RDATA_ERROR_OPEN = 1,
    RDATA_ERROR_SEEK,
    RDATA_ERROR_READ,
    RDATA_ERROR_MALLOC,
    RDATA_ERROR_USER_ABORT,
    RDATA_ERROR_PARSE,
    RDATA_ERROR_WRITE,
    RDATA_ERROR_FACTOR,
    RDATA_ERROR_UNSUPPORTED_COMPRESSION,
    RDATA_ERROR_UNSUPPORTED_CHARSET,
    RDATA_ERROR_CONVERT,
    RDATA_ERROR_CONVERT_BAD_STRING,
    RDATA_ERROR_CONVERT_LONG_STRING,
    RDATA_ERROR_CONVERT_SHORT_STRING,
    RDATA_ERROR_UNSUPPORTED_S_EXPRESSION,
    RDATA_ERROR_UNSUPPORTED_STORAGE_CLASS
} rdata_error_t;

typedef enum rdata_file_format_e {
    RDATA_WORKSPACE,
    RDATA_SINGLE_OBJECT
} rdata_file_format_t;

const char *rdata_error_message(rdata_error_t error_code);

typedef int (*rdata_column_handler)(const char *name, rdata_type_t type,
        void *data, long count, void *ctx);
typedef int (*rdata_table_handler)(const char *name, void *ctx);
typedef int (*rdata_text_value_handler)(
    const char *value, int index, void *ctx
);
typedef int (*rdata_column_name_handler)(
    const char *value, int index, void *ctx
);
typedef void (*rdata_error_handler)(const char *error_message, void *ctx);
typedef int (*rdata_progress_handler)(double progress, void *ctx);

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
typedef __int64 rdata_off_t;
#elif defined _WIN32 || defined __CYGWIN__
typedef _off64_t rdata_off_t;
#elif defined _AIX
typedef off64_t rdata_off_t;
#else
typedef off_t rdata_off_t;
#endif

typedef enum rdata_io_flags_e {
    RDATA_SEEK_SET,
    RDATA_SEEK_CUR,
    RDATA_SEEK_END
} rdata_io_flags_t;

typedef int (*rdata_open_handler)(const char *path, void *io_ctx);
typedef int (*rdata_close_handler)(void *io_ctx);
typedef rdata_off_t (*rdata_seek_handler)(
    rdata_off_t offset, rdata_io_flags_t whence, void *io_ctx
);
typedef ssize_t (*rdata_read_handler)(void *buf, size_t nbyte, void *io_ctx);
typedef rdata_error_t (*rdata_update_handler)(
    long file_size,
    rdata_progress_handler progress_handler,
    void *user_ctx,
    void *io_ctx
);

typedef struct rdata_io_s {
    rdata_open_handler          open;
    rdata_close_handler         close;
    rdata_seek_handler          seek;
    rdata_read_handler          read;
    rdata_update_handler        update;
    void                          *io_ctx;
    int                            external_io;
} rdata_io_t;

typedef struct rdata_parser_s {
    rdata_table_handler         table_handler;
    rdata_column_handler        column_handler;
    rdata_column_name_handler   column_name_handler;
    rdata_column_name_handler   row_name_handler;
    rdata_text_value_handler    text_value_handler;
    rdata_text_value_handler    value_label_handler;
    rdata_column_handler        dim_handler;
    rdata_text_value_handler    dim_name_handler;
    rdata_error_handler         error_handler;
    rdata_io_t                 *io;
} rdata_parser_t;

rdata_parser_t *rdata_parser_init(void);
void rdata_parser_free(rdata_parser_t *parser);

rdata_error_t rdata_set_table_handler(
    rdata_parser_t *parser, rdata_table_handler table_handler
);
rdata_error_t rdata_set_column_handler(
    rdata_parser_t *parser, rdata_column_handler column_handler
);
rdata_error_t rdata_set_column_name_handler(
    rdata_parser_t *parser, rdata_column_name_handler column_name_handler
);
rdata_error_t rdata_set_row_name_handler(
    rdata_parser_t *parser, rdata_column_name_handler row_name_handler
);
rdata_error_t rdata_set_text_value_handler(
    rdata_parser_t *parser, rdata_text_value_handler text_value_handler
);
rdata_error_t rdata_set_value_label_handler(
    rdata_parser_t *parser, rdata_text_value_handler value_label_handler
);
rdata_error_t rdata_set_dim_handler(
    rdata_parser_t *parser, rdata_column_handler dim_handler
);
rdata_error_t rdata_set_dim_name_handler(
    rdata_parser_t *parser, rdata_text_value_handler dim_name_handler
);
rdata_error_t rdata_set_error_handler(
    rdata_parser_t *parser, rdata_error_handler error_handler
);
rdata_error_t rdata_set_open_handler(
    rdata_parser_t *parser, rdata_open_handler open_handler
);
rdata_error_t rdata_set_close_handler(
    rdata_parser_t *parser, rdata_close_handler close_handler
);
rdata_error_t rdata_set_seek_handler(
    rdata_parser_t *parser, rdata_seek_handler seek_handler
);
rdata_error_t rdata_set_read_handler(
    rdata_parser_t *parser, rdata_read_handler read_handler
);
rdata_error_t rdata_set_update_handler(
    rdata_parser_t *parser, rdata_update_handler update_handler);
rdata_error_t rdata_set_io_ctx(rdata_parser_t *parser, void *io_ctx);

/* rdata_parse works on RData and RDS. The table handler will be called once
 * per data frame in RData files, and zero times on RDS files. */

rdata_error_t rdata_parse(
    rdata_parser_t *parser, const char *filename, void *user_ctx
);


// Write API
typedef ssize_t (*rdata_data_writer)(const void *data, size_t len, void *ctx);

typedef struct rdata_column_s {
    rdata_type_t    type;
    int             index;
    char            name[256];
    char            label[1024];

    int32_t         factor_count;
    char          **factor;
} rdata_column_t;

typedef struct rdata_writer_s {
    rdata_file_format_t file_format;
    rdata_data_writer   data_writer;
    size_t              bytes_written;

    rdata_error_handler error_handler;
    void               *user_ctx;

    void               *atom_table;
    int                 bswap;

    rdata_column_t    **columns;
    int32_t             columns_count;
    int32_t             columns_capacity;
} rdata_writer_t;

rdata_writer_t *rdata_writer_init(
    rdata_data_writer write_callback, rdata_file_format_t format
);
void rdata_writer_free(rdata_writer_t *writer);

rdata_column_t *rdata_add_column(
    rdata_writer_t *writer, const char *name, rdata_type_t type
);

rdata_error_t rdata_column_set_label(
    rdata_column_t *column, const char *label
);
rdata_error_t rdata_column_add_factor(
    rdata_column_t *column, const char *factor);

rdata_column_t *rdata_get_column(rdata_writer_t *writer, int32_t j
);

rdata_error_t rdata_begin_file(rdata_writer_t *writer, void *ctx);
rdata_error_t rdata_begin_table(
    rdata_writer_t *writer, const char *variable_name);
rdata_error_t rdata_begin_column(
    rdata_writer_t *writer, rdata_column_t *column, int32_t row_count
);

rdata_error_t rdata_append_real_value(
    rdata_writer_t *writer, double value
);
rdata_error_t rdata_append_int32_value(
    rdata_writer_t *writer, int32_t value
);
rdata_error_t rdata_append_timestamp_value(
    rdata_writer_t *writer, time_t value
);
rdata_error_t rdata_append_date_value(
    rdata_writer_t *writer, struct tm *value
);
rdata_error_t rdata_append_logical_value(
    rdata_writer_t *writer, int value);
rdata_error_t rdata_append_string_value(
    rdata_writer_t *writer, const char *value
);

rdata_error_t rdata_end_column(
    rdata_writer_t *writer, rdata_column_t *column
);
rdata_error_t rdata_end_table(
    rdata_writer_t *writer, int32_t row_count, const char *datalabel
);
rdata_error_t rdata_end_file(rdata_writer_t *writer);

#ifdef __cplusplus
}  // extern c block
#endif

#endif  // PANDAS__LIBS_SRC_LIBRDATA_RDATA_H_
