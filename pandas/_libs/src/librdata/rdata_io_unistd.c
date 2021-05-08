/*
Copyright (c) 2020 Evan Miller
*/

#include <fcntl.h>
#include <stdlib.h>
#if defined _WIN32 || defined __CYGWIN__
    #include <io.h>
#else
    #include <unistd.h>
#endif


#include "rdata.h"
#include "rdata_io_unistd.h"

#if defined _WIN32 || defined __CYGWIN__
#define UNISTD_OPEN_OPTIONS O_RDONLY | O_BINARY
#elif defined _AIX
#define UNISTD_OPEN_OPTIONS O_RDONLY | O_LARGEFILE
#else
#define UNISTD_OPEN_OPTIONS O_RDONLY
#endif

#if defined _WIN32 || defined _AIX
#define lseek lseek
#endif


int rdata_unistd_open_handler(const char *path, void *io_ctx) {
    int fd = open(path, UNISTD_OPEN_OPTIONS);
    ((rdata_unistd_io_ctx_t*) io_ctx)->fd = fd;
    return fd;
}

int rdata_unistd_close_handler(void *io_ctx) {
    int fd = ((rdata_unistd_io_ctx_t*) io_ctx)->fd;
    if (fd != -1)
        return close(fd);
    else
        return 0;
}

rdata_off_t rdata_unistd_seek_handler(
    rdata_off_t offset,
    rdata_io_flags_t whence,
    void *io_ctx
) {
    int flag = 0;
    switch (whence) {
        case RDATA_SEEK_SET:
            flag = SEEK_SET;
            break;
        case RDATA_SEEK_CUR:
            flag = SEEK_CUR;
            break;
        case RDATA_SEEK_END:
            flag = SEEK_END;
            break;
        default:
            return -1;
    }
    int fd = ((rdata_unistd_io_ctx_t*) io_ctx)->fd;
    return lseek(fd, offset, flag);
}

ssize_t rdata_unistd_read_handler(void *buf, size_t nbyte, void *io_ctx) {
    int fd = ((rdata_unistd_io_ctx_t*) io_ctx)->fd;
    ssize_t out = read(fd, buf, nbyte);
    return out;
}

rdata_error_t rdata_unistd_update_handler(long file_size,
        rdata_progress_handler progress_handler, void *user_ctx,
        void *io_ctx) {
    if (!progress_handler)
        return RDATA_OK;

    int fd = ((rdata_unistd_io_ctx_t*) io_ctx)->fd;
    long current_offset = lseek(fd, 0, SEEK_CUR);

    if (current_offset == -1)
        return RDATA_ERROR_SEEK;

    if (progress_handler(1.0 * current_offset / file_size, user_ctx))
        return RDATA_ERROR_USER_ABORT;

    return RDATA_OK;
}

void rdata_unistd_io_init(rdata_parser_t *parser) {
    rdata_set_open_handler(parser, rdata_unistd_open_handler);
    rdata_set_close_handler(parser, rdata_unistd_close_handler);
    rdata_set_seek_handler(parser, rdata_unistd_seek_handler);
    rdata_set_read_handler(parser, rdata_unistd_read_handler);
    rdata_set_update_handler(parser, rdata_unistd_update_handler);

    rdata_unistd_io_ctx_t *io_ctx = calloc(1, sizeof(rdata_unistd_io_ctx_t));
    io_ctx->fd = -1;
    rdata_set_io_ctx(parser, (void*) io_ctx);
}
