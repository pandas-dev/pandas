cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)


cdef extern from "parser/common.h":

    ctypedef enum ParserState:
        START_RECORD
        START_FIELD
        ESCAPED_CHAR
        IN_FIELD
        IN_QUOTED_FIELD
        ESCAPE_IN_QUOTED_FIELD
        QUOTE_IN_QUOTED_FIELD
        EAT_CRNL
        FINISHED

    ctypedef struct table_chunk:
        void **columns
        int ncols

    ctypedef struct parser_t:
        void *source
        char sourcetype   # 'M' for mmap, 'F' for FILE, 'A' for array

        int chunksize  # Number of bytes to prepare for each chunk
        char *data     # pointer to data to be processed
        int datalen    # amount of data available

        # where to write out tokenized data
        char *stream
        int stream_len
        int stream_cap

        # Store words in (potentially ragged) matrix for now, hmm
        char **words
        int *word_starts # where we are in the stream
        int words_len
        int words_cap

        char *pword_start    # pointer to stream start of current field
        int word_start       # position start of current field

        int *line_start      # position in words for start of line
        int *line_fields     # Number of fields in each line
        int lines            # Number of lines observed
        int lines_cap        # Vector capacity

        # Tokenizing stuff
        ParserState state
        int doublequote            # is " represented by ""? */
        char delimiter             # field separator */
        char quotechar             # quote character */
        char escapechar            # escape character */
        int skipinitialspace       # ignore spaces following delimiter? */
        int quoting                # style of quoting to write */

        # hmm =/
        int numeric_field

        char commentchar
        int allow_embedded_newline
        int strict                 # raise exception on bad CSV */

        int error_bad_lines
        int warn_bad_lines

        int infer_types

        # floating point options
        char decimal
        char sci

        # thousands separator (comma, period)
        char thousands

        int header # Boolean: 1: has header, 0: no header

        int skiprows
        int skip_footer

        table_chunk *chunks
        int nchunks

        void **columns
        int ncols

        # PyObject *converters

        #  error handling
        char *error_msg

    parser_t* parser_new()

    int parser_init(parser_t *self)
    void parser_free(parser_t *self)

    void set_parser_default_options(parser_t *self)

    int parser_file_source_init(parser_t *self, FILE* fp)
    int parser_mmap_init(parser_t *self, FILE* fp)
    int parser_array_source_init(parser_t *self, char *bytes, size_t length)

    void debug_print_parser(parser_t *self)

    int tokenize_all_rows(parser_t *self)
    int tokenize_nrows(parser_t *self, size_t nrows)

cdef class TextReader:
    '''

    # source: StringIO or file object

    '''

    cdef:
        parser_t *parser
        object file_handle, should_close

    def __cinit__(self, source, delimiter=',', header=0, memory_map=False,
                  delim_whitespace=False,
                  na_values=None,
                  converters=None,
                  thousands=None):
        self.parser = parser_new()

        self._setup_parser_source(source)

        if delim_whitespace:
            raise NotImplementedError
        else:
            if len(delimiter) > 1:
                raise ValueError('only length-1 separators excluded right now')
            self.parser.delimiter = delimiter

        # TODO: no header vs. header is not the first row
        self.parser.header = header

        self.should_close = False

        self.delimiter = delimiter
        self.delim_whitespace = delim_whitespace

        self.na_values
        self.converters = converters
        self.thousands = thousands

    def __dealloc__(self):
        parser_free(self.parser)

    def __del__(self):
        if self.should_close:
            self.file_handle.close()

    cdef _setup_parser_source(self, source):
        cdef int status

        if isinstance(source, (basestring, file)):
            if isinstance(source, basestring):
                self.file_handle = open(source, 'rb')
                self.should_close = True
                source = self.file_handle

            status = parser_file_source_init(self.parser,
                                             PyFile_AsFile(source))
            if status != 0:
                raise Exception('Initializing from file failed')

        elif hasattr(source, 'read'):
            # e.g., StringIO

            bytes = source.read()

            # TODO: unicode
            if isinstance(bytes, unicode):
                raise ValueError('Only ascii/bytes supported at the moment')

            status = parser_array_source_init(self.parser,
                                              <char*> bytes, len(bytes))
            if status != 0:
                raise Exception('Initializing parser from file-like '
                                'object failed')

    def _parse_table_header(self):
        pass

    def read(self, rows=None):
        """
        rows=None --> read all rows
        """
        cdef int status

        if rows is not None:
            raise NotImplementedError
        else:
            status = tokenize_all_rows(self.parser)


class CParserError(Exception):
    pass


cdef raise_parser_error(object base, parser_t *parser):
    message = '%s. C error: '
    if parser.error_msg != NULL:
        message += parser.error_msg
    else:
        message += 'no error message set'

    raise CParserError(message)
