/*
Copyright (c) 2020 Evan Miller
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "CKHashTable.h"
#include "rdata.h"
#include "rdata_internal.h"

#define R_TAG           0x01
#define R_OBJECT        0x02
#define R_ATTRIBUTES    0x04

#define INITIAL_COLUMNS_CAPACITY    100

#ifdef _WIN32
#define timegm _mkgmtime
#endif

rdata_writer_t *rdata_writer_init(
    rdata_data_writer write_callback,
    rdata_file_format_t format
) {
    rdata_writer_t *writer = calloc(1, sizeof(rdata_writer_t));
    writer->file_format = format;
    writer->bswap = machine_is_little_endian();
    writer->atom_table = ck_hash_table_init(100, 24);
    writer->data_writer = write_callback;

    writer->columns_capacity = INITIAL_COLUMNS_CAPACITY;
    writer->columns = malloc(
        writer->columns_capacity * sizeof(rdata_column_t *));

    return writer;
}

void rdata_writer_free(rdata_writer_t *writer) {
    ck_hash_table_free(writer->atom_table);
    int i, j;
    for (i=0; i < writer->columns_count; i++) {
        rdata_column_t *column = writer->columns[i];
        for (j=0; j < column->factor_count; j++) {
            free(column->factor[j]);
        }
        free(column->factor);
        free(column);
    }
    free(writer->columns);
    free(writer);
}

rdata_column_t *rdata_add_column(
    rdata_writer_t *writer,
    const char *name,
    rdata_type_t type
) {
    if (writer->columns_count == writer->columns_capacity) {
        writer->columns_capacity *= 2;
        writer->columns = realloc(writer->columns,
                writer->columns_capacity * sizeof(rdata_column_t *));
    }
    rdata_column_t *new_column = calloc(1, sizeof(rdata_column_t));

    new_column->index = writer->columns_count++;

    writer->columns[new_column->index] = new_column;

    new_column->type = type;

    if (name) {
        snprintf(new_column->name, sizeof(new_column->name), "%s", name);
    }

    return new_column;
}

rdata_column_t *rdata_get_column(rdata_writer_t *writer, int32_t j) {
    return writer->columns[j];
}

rdata_error_t rdata_column_set_label(
    rdata_column_t *column,
    const char *label
) {
    snprintf(column->label, sizeof(column->label), "%s", label);
    return RDATA_OK;
}

rdata_error_t rdata_column_add_factor(
    rdata_column_t *column,
    const char *factor
) {
    if (column->type != RDATA_TYPE_INT32)
        return RDATA_ERROR_FACTOR;

    char *factor_copy = malloc(strlen(factor)+1);
    strcpy(factor_copy, factor);  // NOLINT

    column->factor_count++;
    column->factor = realloc(
        column->factor,
        sizeof(char *) * column->factor_count);
    column->factor[column->factor_count-1] = factor_copy;

    return RDATA_OK;
}

static rdata_error_t rdata_write_bytes(
    rdata_writer_t *writer,
    const void *data, size_t len
) {
    size_t bytes_written = writer->data_writer(data, len, writer->user_ctx);
    if (bytes_written < len) {
        return RDATA_ERROR_WRITE;
    }
    writer->bytes_written += bytes_written;
    return RDATA_OK;
}

static rdata_error_t rdata_write_integer(
    rdata_writer_t *writer,
    int32_t val
) {
    if (writer->bswap) {
        val = byteswap4(val);
    }
    return rdata_write_bytes(writer, &val, sizeof(val));
}

static rdata_error_t rdata_write_double(rdata_writer_t *writer, double val) {
    if (writer->bswap) {
        val = byteswap_double(val);
    }
    return rdata_write_bytes(writer, &val, sizeof(val));
}

static rdata_error_t rdata_write_header(
    rdata_writer_t *writer,
    int type,
    int flags
) {
    rdata_sexptype_header_t header;
    memset(&header, 0, sizeof(header));

    header.type = type;
    header.object = !!(flags & R_OBJECT);
    header.tag = !!(flags & R_TAG);
    header.attributes = !!(flags & R_ATTRIBUTES);

    uint32_t sexp_int;

    memcpy(&sexp_int, &header, sizeof(header));

    return rdata_write_integer(writer, sexp_int);
}

static rdata_error_t rdata_write_string(
    rdata_writer_t *writer,
    const char *string
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_header(writer, RDATA_SEXPTYPE_CHARACTER_STRING, 0);
    if (retval != RDATA_OK)
        goto cleanup;

    ssize_t len = string ? (ssize_t)strlen(string) : -1;

    retval = rdata_write_integer(writer, len);
    if (retval != RDATA_OK)
        goto cleanup;

    if (len > 0)
        return rdata_write_bytes(writer, string, len);

cleanup:
    return retval;
}

static rdata_error_t rdata_write_pairlist_key(
    rdata_writer_t *writer,
    const char *key
) {
    rdata_error_t retval = RDATA_OK;
    ck_hash_table_t *atom_table = (ck_hash_table_t *)writer->atom_table;
    uint64_t ref = (uint64_t)ck_str_hash_lookup(key, atom_table);
    if (ref == 0) {
        ck_str_hash_insert(key, (void *)(atom_table->count + 1), atom_table);

        retval = rdata_write_integer(writer, 1);
        if (retval != RDATA_OK)
            goto cleanup;

        retval = rdata_write_string(writer, key);
    } else {
        retval = rdata_write_integer(writer, (ref << 8) | 0xFF);
    }

cleanup:
    return retval;
}

static rdata_error_t rdata_write_pairlist_header(
    rdata_writer_t *writer,
    const char *key
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_header(writer, RDATA_SEXPTYPE_PAIRLIST, R_TAG);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_pairlist_key(writer, key);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_write_attributed_vector_header(
    rdata_writer_t *writer, int type,
    int32_t size
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_header(writer, type, R_OBJECT | R_ATTRIBUTES);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_integer(writer, size);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_write_simple_vector_header(
    rdata_writer_t *writer,
    int type,
    int32_t size
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_header(writer, type, 0);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_integer(writer, size);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_write_class_pairlist(
    rdata_writer_t *writer,
    const char *class
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_pairlist_header(writer, "class");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_CHARACTER_VECTOR,
        1);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_string(writer, class);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

rdata_error_t rdata_begin_file(
    rdata_writer_t *writer,
    void *user_ctx
) {
    rdata_error_t retval = RDATA_OK;

    writer->user_ctx = user_ctx;

    if (writer->file_format == RDATA_WORKSPACE) {
        retval = rdata_write_bytes(writer, "RDX2\n", 5);
        if (retval != RDATA_OK)
            goto cleanup;
    }

    rdata_v2_header_t v2_header;
    memcpy(v2_header.header, "X\n", sizeof("X\n")-1);
    v2_header.format_version = 2;
    v2_header.reader_version = 131840;
    v2_header.writer_version = 131840;

    if (writer->bswap) {
        v2_header.format_version = byteswap4(v2_header.format_version);
        v2_header.reader_version = byteswap4(v2_header.reader_version);
        v2_header.writer_version = byteswap4(v2_header.writer_version);
    }

    retval = rdata_write_bytes(writer, &v2_header, sizeof(v2_header));
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

rdata_error_t rdata_begin_table(
    rdata_writer_t *writer,
    const char *variable_name
) {
    rdata_error_t retval = RDATA_OK;

    if (writer->file_format == RDATA_WORKSPACE) {
        retval = rdata_write_pairlist_header(writer, variable_name);
        if (retval != RDATA_OK)
            goto cleanup;
    }

    retval = rdata_write_attributed_vector_header(
         writer,
         RDATA_SEXPTYPE_GENERIC_VECTOR,
         writer->columns_count);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_begin_factor_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_attributed_vector_header(
        writer,
        RDATA_SEXPTYPE_INTEGER_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_factor_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    int i;

    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_pairlist_header(writer, "levels");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_simple_vector_header(writer,
            RDATA_SEXPTYPE_CHARACTER_VECTOR, column->factor_count);
    if (retval != RDATA_OK)
        goto cleanup;

    for (i=0; i < column->factor_count; i++) {
        retval = rdata_write_string(writer, column->factor[i]);
        if (retval != RDATA_OK)
            goto cleanup;
    }

    retval = rdata_write_class_pairlist(writer, "factor");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_header(writer, RDATA_PSEUDO_SXP_NIL, 0);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_begin_real_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_REAL_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_real_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    return RDATA_OK;
}

static rdata_error_t rdata_begin_timestamp_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_attributed_vector_header(
        writer,
        RDATA_SEXPTYPE_REAL_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_timestamp_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_class_pairlist(writer, "POSIXct");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_header(writer, RDATA_PSEUDO_SXP_NIL, 0);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_begin_date_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_attributed_vector_header(
        writer,
        RDATA_SEXPTYPE_REAL_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_date_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_class_pairlist(writer, "Date");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_header(writer, RDATA_PSEUDO_SXP_NIL, 0);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

static rdata_error_t rdata_begin_integer_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_INTEGER_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_integer_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    return RDATA_OK;
}

static rdata_error_t rdata_begin_logical_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_LOGICAL_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_logical_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    return RDATA_OK;
}

static rdata_error_t rdata_begin_string_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    return rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_CHARACTER_VECTOR,
        row_count);
}

static rdata_error_t rdata_end_string_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    return RDATA_OK;
}

rdata_error_t rdata_begin_column(
    rdata_writer_t *writer,
    rdata_column_t *column,
    int32_t row_count
) {
    rdata_type_t type = column->type;

    if (type == RDATA_TYPE_INT32) {
        if (column->factor_count)
            return rdata_begin_factor_column(writer, column, row_count);
        return rdata_begin_integer_column(writer, column, row_count);
    }
    if (type == RDATA_TYPE_REAL)
        return rdata_begin_real_column(writer, column, row_count);
    if (type == RDATA_TYPE_TIMESTAMP)
        return rdata_begin_timestamp_column(writer, column, row_count);
    if (type == RDATA_TYPE_DATE)
        return rdata_begin_date_column(writer, column, row_count);
    if (type == RDATA_TYPE_LOGICAL)
        return rdata_begin_logical_column(writer, column, row_count);
    if (type == RDATA_TYPE_STRING)
        return rdata_begin_string_column(writer, column, row_count);

    return RDATA_OK;
}

rdata_error_t rdata_append_real_value(
    rdata_writer_t *writer,
    double value
) {
    return rdata_write_double(writer, value);
}

rdata_error_t rdata_append_int32_value(
    rdata_writer_t *writer,
    int32_t value
) {
    return rdata_write_integer(writer, value);
}

rdata_error_t rdata_append_timestamp_value(
    rdata_writer_t *writer,
    time_t value
) {
    return rdata_write_double(writer, value);
}

rdata_error_t rdata_append_date_value(
    rdata_writer_t *writer,
    struct tm *value
) {
    return rdata_write_double(writer, timegm(value) / 86400);
}

rdata_error_t rdata_append_logical_value(
    rdata_writer_t *writer,
    int value
) {
    if (value < 0)
        return rdata_write_integer(writer, INT32_MIN);

    return rdata_write_integer(writer, (value > 0));
}

rdata_error_t rdata_append_string_value(
    rdata_writer_t *writer,
    const char *value
) {
    return rdata_write_string(writer, value);
}

rdata_error_t rdata_end_column(
    rdata_writer_t *writer,
    rdata_column_t *column
) {
    rdata_type_t type = column->type;

    if (type == RDATA_TYPE_INT32) {
        if (column->factor_count)
            return rdata_end_factor_column(writer, column);
        return rdata_end_integer_column(writer, column);
    }
    if (type == RDATA_TYPE_REAL)
        return rdata_end_real_column(writer, column);
    if (type == RDATA_TYPE_TIMESTAMP)
        return rdata_end_timestamp_column(writer, column);
    if (type == RDATA_TYPE_DATE)
        return rdata_end_date_column(writer, column);
    if (type == RDATA_TYPE_LOGICAL)
        return rdata_end_logical_column(writer, column);
    if (type == RDATA_TYPE_STRING)
        return rdata_end_string_column(writer, column);

    return RDATA_OK;
}

rdata_error_t rdata_end_table(
    rdata_writer_t *writer,
    int32_t row_count,
    const char *datalabel
) {
    int i;
    rdata_error_t retval = RDATA_OK;

    retval = rdata_write_pairlist_header(writer, "datalabel");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_simple_vector_header(
        writer,
        RDATA_SEXPTYPE_CHARACTER_VECTOR,
        1);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_string(writer, datalabel);
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_pairlist_header(writer, "names");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_simple_vector_header(writer,
            RDATA_SEXPTYPE_CHARACTER_VECTOR, writer->columns_count);
    if (retval != RDATA_OK)
        goto cleanup;

    for (i=0; i < writer->columns_count; i++) {
        retval = rdata_write_string(writer, writer->columns[i]->name);
        if (retval != RDATA_OK)
            goto cleanup;
    }

    retval = rdata_write_pairlist_header(writer, "var.labels");
    if (retval != RDATA_OK)
        goto cleanup;

    retval = rdata_write_simple_vector_header(writer,
            RDATA_SEXPTYPE_CHARACTER_VECTOR, writer->columns_count);
    if (retval != RDATA_OK)
        goto cleanup;

    for (i=0; i < writer->columns_count; i++) {
        retval = rdata_write_string(writer, writer->columns[i]->label);
        if (retval != RDATA_OK)
            goto cleanup;
    }

    retval = rdata_write_class_pairlist(writer, "data.frame");
    if (retval != RDATA_OK)
        goto cleanup;

    if (row_count > 0) {
        retval = rdata_write_pairlist_header(writer, "row.names");
        if (retval != RDATA_OK)
            goto cleanup;

        retval = rdata_write_simple_vector_header(writer,
                RDATA_SEXPTYPE_CHARACTER_VECTOR, row_count);
        if (retval != RDATA_OK)
            goto cleanup;

        char buf[128];
        for (i=0; i < row_count; i++) {
            snprintf(buf, sizeof(buf), "%d", i+1);
            retval = rdata_write_string(writer, buf);
            if (retval != RDATA_OK)
                goto cleanup;
        }
    }

    retval = rdata_write_header(writer, RDATA_PSEUDO_SXP_NIL, 0);
    if (retval != RDATA_OK)
        goto cleanup;

cleanup:
    return retval;
}

rdata_error_t rdata_end_file(rdata_writer_t *writer) {
    if (writer->file_format == RDATA_WORKSPACE)
        return rdata_write_header(writer, RDATA_PSEUDO_SXP_NIL, 0);

    return RDATA_OK;
}
