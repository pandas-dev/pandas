# cython: c_string_type=str, c_string_encoding=utf8, language_level=3

cdef int handle_table(const char *name, void *ctx) except *:
    """
    Retrieves original R object name.

    Called once per data frame in RData files,
    and zero times on RDS files.
    """
    lbr = <LibrdataReader>ctx

    lbr.colidx = 0
    lbr.rows = 0
    lbr.rlevels = {}
    lbr.rtext = {}
    lbr.is_factor = False
    lbr.rownames = {}
    lbr.colnames = {}
    lbr.dims = 0
    lbr.dim_str = {}

    if name != NULL:
        lbr.tblname = name

    if "r_dataframe" in lbr.rvalues.keys():
        lbr.rvalues[lbr.tblname] = lbr.rvalues.pop("r_dataframe")
    else:
        lbr.rvalues[lbr.tblname] = {
            "data": {},
            "dtypes": {},
            "colnames": None,
            "rownames": None
        }
    return 0     # non-zero to abort processing


cdef int handle_column(
    const char *name,
    rdata_type_t dtype,
    void *data,
    long count,
    void *ctx
) except *:
    """
    Parses each non-string column in data frame.

    Called once for all columns with the following caveats:
    * `name` is NULL for some columns (see handle_column_name below)
    * `data` is NULL for text columns (see handle_text_value below)
    Special conditon for matrices with dims attribute.
    """
    lbr = <LibrdataReader>ctx

    lbr.rows = count
    cdef int *rints = <int*>data
    cdef double *rdoubles = <double*>data

    if dtype in [
        rdata_type_t.RDATA_TYPE_REAL,
        rdata_type_t.RDATA_TYPE_DATE,
        rdata_type_t.RDATA_TYPE_TIMESTAMP
    ]:
        lbr.rvalues[lbr.tblname]["dtypes"][lbr.colidx] = lbr.rtypes[dtype]
        lbr.rvalues[lbr.tblname]["data"][lbr.colidx] = {
            i: rdoubles[i] for i in range(count)
        }
        lbr.colidx += 1

    elif dtype in [
        rdata_type_t.RDATA_TYPE_INT32,
        rdata_type_t.RDATA_TYPE_LOGICAL
    ]:
        if lbr.is_factor:
            lbr.rvalues[lbr.tblname]["dtypes"][lbr.colidx] = "factor"
            lbr.rvalues[lbr.tblname]["data"][lbr.colidx] = {
                i: float('nan') if rints[i] < 0 else lbr.rlevels[rints[i]-1]
                for i in range(count)
            }
            lbr.is_factor = False
        else:
            lbr.rvalues[lbr.tblname]["dtypes"][lbr.colidx] = lbr.rtypes[dtype]
            lbr.rvalues[lbr.tblname]["data"][lbr.colidx] = {
                i: rints[i] for i in range(count)
            }
        lbr.colidx += 1

    if lbr.dims > 0:
        lbr.tblname = "r_matrix"
        lbr.rvalues[lbr.tblname] = lbr.rvalues.pop("r_dataframe")
        dim_data = list(lbr.rvalues[lbr.tblname]["data"][0].values())

        n = 0
        rows, cols = lbr.dim_str.values()
        for col in range(cols):
            lbr.rvalues[lbr.tblname]["dtypes"][col] = lbr.rtypes[dtype]
            lbr.rvalues[lbr.tblname]["data"][col] = {
                i: d for i, d in enumerate(dim_data[n:n+rows])
            }
            n += rows

    return 0

cdef int handle_text_value(const char *value, int index, void *ctx) except *:
    """
    Parses string data.

    Called once per row for a text column.
    """
    lbr = <LibrdataReader>ctx

    lbr.rtext[index] = value if value != NULL else None

    if index == (lbr.rows - 1):
        lbr.rvalues[lbr.tblname]["dtypes"][lbr.colidx] = "str"
        lbr.rvalues[lbr.tblname]["data"][lbr.colidx] = lbr.rtext
        lbr.colidx += 1
        lbr.rtext = {}

    return 0

cdef int handle_value_label(const char *value, int index, void *ctx) except *:
    """
    Parses factor levels.

    Called for factor variables, once for each level
    """
    lbr = <LibrdataReader>ctx

    lbr.is_factor = True
    lbr.rlevels[index] = value

    return 0

cdef int handle_dim(
    const char *name,
    rdata_type_t dtype,
    void *data,
    long count,
    void *ctx
) except *:
    """
    Parses meta data on non-dataframe objects

    Called once for objects with R dims (matrices, arrays, etc.)).
    Special condition for character matrices.
    """
    lbr = <LibrdataReader>ctx

    cdef int *rdims = <int*>data

    lbr.dims = count
    lbr.dim_str = {i: rdims[i] for i in range(count)}

    if lbr.rvalues[lbr.tblname]["dtypes"] == {0: "str"}:
        dim_data = list(lbr.rvalues[lbr.tblname]["data"][0].values())

        n = 0
        rows, cols = lbr.dim_str.values()

        for col in range(cols):
            lbr.rvalues[lbr.tblname]["dtypes"][col] = "str"
            lbr.rvalues[lbr.tblname]["data"][col] = dim_data[n:n+rows]
            n += rows

    return 0

cdef int handle_column_name(const char *name, int index, void *ctx) except *:
    """
    Retrieves column names of data frame

    Returns only non-NULL column names after parsing data.
    """
    lbr = <LibrdataReader>ctx

    lbr.colnames[index] = name
    lbr.rvalues[lbr.tblname]["colnames"] = lbr.colnames

    return 0

cdef int handle_row_name(const char *name, int index, void *ctx) except *:
    """
    Retrieves row names of data frame

    Returns only non-NULL row names appear after parsing data.
    """
    lbr = <LibrdataReader>ctx

    lbr.rownames[index] = name
    lbr.rvalues[lbr.tblname]["rownames"] = lbr.rownames

    return 0

cdef int handle_dim_name(const char *name, int index, void *ctx) except *:
    """
    Retrieves dim names of matrices or arrays

    Returns only non-NULL dim names appear after parsing data.
    """

    lbr = <LibrdataReader>ctx

    if (index < lbr.dim_str[0]) and lbr.rownames.get(index) is None:
        lbr.rownames[index] = name if name != NULL else str(index)
    else:
        lbr.rvalues[lbr.tblname]["rownames"] = lbr.rownames

    if index < lbr.dim_str[1]:
        lbr.colnames[index] = name if name != NULL else str(index)
    else:
        lbr.rvalues[lbr.tblname]["colnames"] = lbr.colnames

    return 0


class LibrdataParserError(Exception):
    """
    Base error class to capture exceptions in librdata parsing.
    """
    pass


cdef class LibrdataReader:
    """
    Base class to read RData files.

    Class interfaces with librdata C library to builds dictionaries
    of each data frame including data content and meta (dtypes, colnames,
    and rownames). Callbacks above are used in ``rdata_`` method attributes.
    """
    cdef rdata_parser_t *rparser
    cdef public:
        int colidx
        int rows
        dict rlevels
        dict rtext
        bint is_factor
        dict rownames
        dict colnames
        dict rtypes
        str tblname
        dict rvalues
        int dims
        dict dim_str

    cpdef read_rdata(self, rfile):
        self.rparser = rdata_parser_init()

        self.colidx = 0
        self.rows = 0
        self.rlevels = {}
        self.rtext = {}
        self.is_factor = False
        self.rownames = {}
        self.colnames = {}
        self.dims = 0
        self.dim_str = {}
        self.rtypes = {
            rdata_type_t.RDATA_TYPE_LOGICAL: "bool",
            rdata_type_t.RDATA_TYPE_INT32: "int",
            rdata_type_t.RDATA_TYPE_REAL: "float",
            rdata_type_t.RDATA_TYPE_DATE: "date",
            rdata_type_t.RDATA_TYPE_TIMESTAMP: "datetime",
            rdata_type_t.RDATA_TYPE_STRING: "str"
        }
        self.tblname = "r_dataframe"
        self.rvalues = {
            "r_dataframe": {
                "data": {},
                "dtypes": {},
                "colnames": None,
                "rownames": None
            }
        }

        err = RDATA_OK
        while err == RDATA_OK:
            err = rdata_set_table_handler(self.rparser, handle_table)
            err = rdata_set_dim_handler(self.rparser, handle_dim)
            err = rdata_set_column_handler(self.rparser, handle_column)
            err = rdata_set_text_value_handler(self.rparser, handle_text_value)
            err = rdata_set_value_label_handler(self.rparser, handle_value_label)
            err = rdata_set_column_name_handler(self.rparser, handle_column_name)
            err = rdata_set_row_name_handler(self.rparser, handle_row_name)
            err = rdata_set_dim_name_handler(self.rparser, handle_dim_name)

            err = rdata_parse(self.rparser, rfile, <void *>self)
            rdata_parser_free(self.rparser)
            break

        if err != RDATA_OK:
            msg = rdata_error_message(err)
            raise LibrdataParserError(msg)

        return self.rvalues


class LibrdataWriterError(Exception):
    """
    Base error class to capture exceptions in librdata writing.
    """
    pass


cdef ssize_t write_data(const void *bytes, size_t len, void *ctx):
    cdef int fd = (<int*>ctx)[0]
    return write(fd, bytes, len)

cdef class LibrdataWriter():
    """
    Base class to write RData files.

    Class interfaces with librdata C library to iterate through dictionaries
    of each DataFrame column according to correspoinding dtype.
    Single callback above is usedd in exposed `init`` method.
    """
    cdef:
        int fd
        int row_count
        dict rdict
        dict rformats
        dict rtypes
        str tbl_name
        rdata_writer_t *writer
        rdata_column_t *py_col

    cdef write_col_data(self, i, kdata, vdata, ktype, vtype):
        py_col = rdata_get_column(self.writer, i)
        rdata_begin_column(self.writer, py_col, self.row_count)

        if vtype == "bool":
            for k, v in vdata.items():
                rdata_append_logical_value(self.writer, v)

        if vtype.startswith(("int", "uint")):
            for k, v in vdata.items():
                rdata_append_int32_value(self.writer, v)

        if vtype.startswith("float"):
            for k, v in vdata.items():
                rdata_append_real_value(self.writer, v)

        if vtype.startswith("datetime64"):
            for k, v in vdata.items():
                rdata_append_timestamp_value(self.writer, v)

        if vtype == "object":
            for k, v in vdata.items():
                if v == v:
                    rdata_append_string_value(self.writer, v)
                else:
                    rdata_append_string_value(self.writer, NULL)

        rdata_end_column(self.writer, py_col)

    cpdef write_rdata(self, rfile, rdict, rformat, tbl_name=None):

        self.rdict = rdict
        self.tbl_name = tbl_name
        self.row_count = len(next(iter(rdict["data"].items()))[1])

        self.rformats = {
            "rdata": RDATA_WORKSPACE,
            "rda": RDATA_WORKSPACE,
            "rds": RDATA_SINGLE_OBJECT
        }

        self.rtypes = {
            "bool": RDATA_TYPE_LOGICAL,
            "int8": RDATA_TYPE_INT32,
            "int16": RDATA_TYPE_INT32,
            "int32": RDATA_TYPE_INT32,
            "int64": RDATA_TYPE_INT32,
            "uint8": RDATA_TYPE_INT32,
            "uint16": RDATA_TYPE_INT32,
            "uint32": RDATA_TYPE_INT32,
            "uint64": RDATA_TYPE_INT32,
            "float8": RDATA_TYPE_REAL,
            "float16": RDATA_TYPE_REAL,
            "float32": RDATA_TYPE_REAL,
            "float64": RDATA_TYPE_REAL,
            "datetime64[ns]": RDATA_TYPE_TIMESTAMP,
            "object": RDATA_TYPE_STRING
        }

        self.fd = open(rfile, O_CREAT | O_WRONLY, 0644);
        self.writer = rdata_writer_init(write_data, self.rformats[rformat])

        for k, v in self.rdict["dtypes"].items():
            rdata_add_column(self.writer, k, self.rtypes[v])

        rdata_begin_file(self.writer, &self.fd)
        rdata_begin_table(self.writer, self.tbl_name)

        try:
            for n, ((kd, vd), (kt, vt)) in enumerate(
                zip(
                    self.rdict["data"].items(),
                    self.rdict["dtypes"].items()
                )
            ):
                self.write_col_data(n, kd, vd, kt, vt)

        except (TypeError, ValueError, UnicodeDecodeError) as e:
            raise LibrdataWriterError(
                "DataFrame contains one more invalid types or data values. "
                "that does not conform to R data types."
            )

        rdata_end_table(self.writer, self.row_count, "pandas_dataframe")
        rdata_end_file(self.writer)

        close(self.fd)
        rdata_writer_free(self.writer)
