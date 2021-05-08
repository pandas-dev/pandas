/*
Copyright (c) 2020 Evan Miller
*/

#include <stdlib.h>
#include "rdata.h"
#include "rdata_io_unistd.h"

rdata_parser_t *rdata_parser_init() {
    rdata_parser_t *parser = calloc(1, sizeof(rdata_parser_t));
    parser->io = calloc(1, sizeof(rdata_io_t));
    rdata_unistd_io_init(parser);
    return parser;
}

void rdata_parser_free(rdata_parser_t *parser) {
    if (parser) {
        if (parser->io)
            free(parser->io);
        free(parser);
    }
}

rdata_error_t rdata_set_table_handler(
    rdata_parser_t *parser,
    rdata_table_handler table_handler
) {
    parser->table_handler = table_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_column_handler(
    rdata_parser_t *parser,
    rdata_column_handler column_handler
) {
    parser->column_handler = column_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_column_name_handler(
    rdata_parser_t *parser,
    rdata_column_name_handler column_name_handler
) {
    parser->column_name_handler = column_name_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_row_name_handler(
    rdata_parser_t *parser,
    rdata_column_name_handler row_name_handler
) {
    parser->row_name_handler = row_name_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_text_value_handler(
    rdata_parser_t *parser,
    rdata_text_value_handler text_value_handler
) {
    parser->text_value_handler = text_value_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_value_label_handler(
    rdata_parser_t *parser,
    rdata_text_value_handler value_label_handler
) {
    parser->value_label_handler = value_label_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_dim_handler(
    rdata_parser_t *parser,
    rdata_column_handler dim_handler
) {
    parser->dim_handler = dim_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_dim_name_handler(
    rdata_parser_t *parser,
    rdata_text_value_handler dim_name_handler
) {
    parser->dim_name_handler = dim_name_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_error_handler(
    rdata_parser_t *parser,
    rdata_error_handler error_handler
) {
    parser->error_handler = error_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_open_handler(
    rdata_parser_t *parser,
    rdata_open_handler open_handler
) {
    parser->io->open = open_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_close_handler(
    rdata_parser_t *parser,
    rdata_close_handler close_handler
) {
    parser->io->close = close_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_seek_handler(
    rdata_parser_t *parser,
    rdata_seek_handler seek_handler
) {
    parser->io->seek = seek_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_read_handler(
    rdata_parser_t *parser,
    rdata_read_handler read_handler
) {
    parser->io->read = read_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_update_handler(
    rdata_parser_t *parser,
    rdata_update_handler update_handler
) {
    parser->io->update = update_handler;
    return RDATA_OK;
}

rdata_error_t rdata_set_io_ctx(
    rdata_parser_t *parser,
    void *io_ctx
) {
    if (!parser->io->external_io)
        free(parser->io->io_ctx);

    parser->io->io_ctx = io_ctx;
    parser->io->external_io = 1;

    return RDATA_OK;
}
