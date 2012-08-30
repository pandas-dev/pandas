
#define READ_ERROR_OUT_OF_MEMORY   1

int count_rows(FILE *f, char delimiter, char quote, char comment,
               int allow_embedded_newline);

int count_fields(FILE *f, char delimiter, char quote, char comment,
                 int allow_embedded_newline);


void *read_rows(FILE *f, int *nrows, char *fmt,
                char delimiter, char quote, char comment,
                char sci, char decimal,
                int allow_embedded_newline,
                char *datetime_fmt,
                int tz_offset,
                int *usecols, int num_usecols,
                int skiprows,
                void *data_array,
                int *p_error_type, int *p_error_lineno);
