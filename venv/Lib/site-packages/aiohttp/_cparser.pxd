# This file contains just the definitions needed in our Cython code.

from libc.stdint cimport uint8_t, uint16_t, uint64_t


cdef extern from "llhttp.h":

    struct llhttp__internal_s:
        void* data
        uint64_t content_length
        uint8_t type
        uint8_t method
        uint8_t http_major
        uint8_t http_minor
        uint8_t upgrade
        uint16_t flags
        uint16_t status_code

    ctypedef llhttp__internal_s llhttp__internal_t
    ctypedef llhttp__internal_t llhttp_t

    ctypedef int (*llhttp_data_cb)(llhttp_t*, const char *at, size_t length) except -1
    ctypedef int (*llhttp_cb)(llhttp_t*) except -1

    struct llhttp_settings_s:
        llhttp_cb      on_message_begin
        llhttp_data_cb on_url
        llhttp_data_cb on_status
        llhttp_data_cb on_header_field
        llhttp_data_cb on_header_value
        llhttp_cb      on_headers_complete
        llhttp_data_cb on_body
        llhttp_cb      on_message_complete
        llhttp_cb      on_chunk_header
        llhttp_cb      on_chunk_complete

    ctypedef llhttp_settings_s llhttp_settings_t

    enum llhttp_errno:
        HPE_OK,
        HPE_INVALID_METHOD,
        HPE_INVALID_URL,
        HPE_INVALID_CONSTANT,
        HPE_INVALID_VERSION,
        HPE_INVALID_HEADER_TOKEN,
        HPE_INVALID_CONTENT_LENGTH,
        HPE_INVALID_CHUNK_SIZE,
        HPE_INVALID_STATUS,
        HPE_INVALID_EOF_STATE,
        HPE_INVALID_TRANSFER_ENCODING,
        HPE_CB_MESSAGE_BEGIN,
        HPE_CB_HEADERS_COMPLETE,
        HPE_CB_MESSAGE_COMPLETE,
        HPE_CB_CHUNK_HEADER,
        HPE_CB_CHUNK_COMPLETE,
        HPE_PAUSED,
        HPE_PAUSED_UPGRADE

    ctypedef llhttp_errno llhttp_errno_t

    enum llhttp_flags:
        F_CHUNKED,
        F_CONTENT_LENGTH

    enum llhttp_type:
        HTTP_REQUEST,
        HTTP_RESPONSE

    enum llhttp_method:
        HTTP_CONNECT

    void llhttp_settings_init(llhttp_settings_t* settings)
    void llhttp_init(llhttp_t* parser, llhttp_type type,
                 const llhttp_settings_t* settings)

    llhttp_errno_t llhttp_execute(llhttp_t* parser, const char* data, size_t len)

    int llhttp_should_keep_alive(const llhttp_t* parser)

    void llhttp_resume(llhttp_t* parser)
    void llhttp_resume_after_upgrade(llhttp_t* parser)

    llhttp_errno_t llhttp_get_errno(const llhttp_t* parser)
    const char* llhttp_get_error_reason(const llhttp_t* parser)
    const char* llhttp_get_error_pos(const llhttp_t* parser)

    void llhttp_set_lenient_headers(llhttp_t* parser, int enabled)
    void llhttp_set_lenient_optional_cr_before_lf(llhttp_t* parser, int enabled)
    void llhttp_set_lenient_spaces_after_chunk_size(llhttp_t* parser, int enabled)
