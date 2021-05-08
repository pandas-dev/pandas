/*
Copyright (c) 2020 Evan Miller
*/

#ifndef PANDAS_IO_RDATA_LIBRDATA_RDATA_IO_UNISTD_H_
#define PANDAS_IO_RDATA_LIBRDATA_RDATA_IO_UNISTD_H_

typedef struct rdata_unistd_io_ctx_s {
    int               fd;
} rdata_unistd_io_ctx_t;

int rdata_unistd_open_handler(const char *path, void *io_ctx);
int rdata_unistd_close_handler(void *io_ctx);
rdata_off_t rdata_unistd_seek_handler(
    rdata_off_t offset, rdata_io_flags_t whence, void *io_ctx
);
ssize_t rdata_unistd_read_handler(void *buf, size_t nbytes, void *io_ctx);
rdata_error_t rdata_unistd_update_handler(
    long file_size,
    rdata_progress_handler progress_handler,
    void *user_ctx,
    void *io_ctx
);
void rdata_unistd_io_init(rdata_parser_t *parser);

#endif  // PANDAS_IO_RDATA_LIBRDATA_RDATA_IO_UNISTD_H_
