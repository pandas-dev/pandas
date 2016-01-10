import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t

# rle_decompress decompresses data using a Run Length Encoding
# algorithm.  It is partially documented here:
#
# https://cran.r-project.org/web/packages/sas7bdat/vignettes/sas7bdat.pdf
def _rle_decompress(int result_length, np.ndarray[uint8_t, ndim=1] inbuff):

    cdef uint8_t control_byte
    cdef uint8_t [:] result = np.zeros(result_length, np.uint8)

    cdef int rpos = 0
    cdef int ipos = 0
    cdef int i
    cdef int nbytes
    cdef uint8_t x
    cdef length = len(inbuff)

    while ipos < length:
        control_byte = inbuff[ipos] & 0xF0
        end_of_first_byte = int(inbuff[ipos] & 0x0F)
        ipos += 1

        if control_byte == 0x00:
            if end_of_first_byte != 0:
                print("Unexpected non-zero end_of_first_byte")
            nbytes = int(inbuff[ipos]) + 64
            ipos += 1
            for i in range(nbytes):
                result[rpos] = inbuff[ipos]
                rpos += 1
                ipos += 1
        elif control_byte == 0x40:
            # not documented
            nbytes = end_of_first_byte * 16
            nbytes += int(inbuff[ipos])
            ipos += 1
            for i in range(nbytes):
                result[rpos] = inbuff[ipos]
                rpos += 1
            ipos += 1
        elif control_byte == 0x60:
            nbytes = end_of_first_byte*256 + int(inbuff[ipos]) + 17
            ipos += 1
            for i in range(nbytes):
                result[rpos] = 0x20
                rpos += 1
        elif control_byte == 0x70:
            nbytes = end_of_first_byte*256 + int(inbuff[ipos]) + 17
            ipos += 1
            for i in range(nbytes):
                result[rpos] = 0x00
                rpos += 1
        elif control_byte == 0x80:
            nbytes = end_of_first_byte + 1
            for i in range(nbytes):
                result[rpos] = inbuff[ipos + i]
                rpos += 1
            ipos += nbytes
        elif control_byte == 0x90:
            nbytes = end_of_first_byte + 17
            for i in range(nbytes):
                result[rpos] = inbuff[ipos + i]
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xA0:
            nbytes = end_of_first_byte + 33
            for i in range(nbytes):
                result[rpos] = inbuff[ipos + i]
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xB0:
            nbytes = end_of_first_byte + 49
            for i in range(nbytes):
                result[rpos] = inbuff[ipos + i]
                rpos += 1
            ipos += nbytes
        elif control_byte == 0xC0:
            nbytes = end_of_first_byte + 3
            x = inbuff[ipos]
            ipos += 1
            for i in range(nbytes):
                result[rpos] = x
                rpos += 1
        elif control_byte == 0xD0:
            nbytes = end_of_first_byte + 2
            for i in range(nbytes):
                result[rpos] = 0x40
                rpos += 1
        elif control_byte == 0xE0:
            nbytes = end_of_first_byte + 2
            for i in range(nbytes):
                result[rpos] = 0x20
                rpos += 1
        elif control_byte == 0xF0:
            nbytes = end_of_first_byte + 2
            for i in range(nbytes):
                result[rpos] = 0x00
                rpos += 1
        else:
            raise ValueError("unknown control byte: %v", control_byte)

    if len(result) != result_length:
        print("RLE: %v != %v\n", (len(result), result_length))

    return np.asarray(result).tostring()


# rdc_decompress decompresses data using the Ross Data Compression algorithm:
#
#   http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
def _rdc_decompress(int result_length, np.ndarray[uint8_t, ndim=1] inbuff):

    cdef uint8_t cmd
    cdef uint16_t ctrl_bits
    cdef uint16_t ctrl_mask = 0
    cdef uint16_t ofs
    cdef uint16_t cnt
    cdef int ipos = 0
    cdef int rpos = 0
    cdef int k

    cdef uint8_t [:] outbuff = np.zeros(result_length, dtype=np.uint8)

    ii = -1

    while ipos < len(inbuff):
        ii += 1
        ctrl_mask = ctrl_mask >> 1
        if ctrl_mask == 0:
            ctrl_bits = (<uint16_t>inbuff[ipos] << 8) + <uint16_t>inbuff[ipos + 1]
            ipos += 2
            ctrl_mask = 0x8000

        if ctrl_bits & ctrl_mask == 0:
            outbuff[rpos] = inbuff[ipos]
            ipos += 1
            rpos += 1
            continue

        cmd = (inbuff[ipos] >> 4) & 0x0F
        cnt = <uint16_t>(inbuff[ipos] & 0x0F)
        ipos += 1

        # short RLE
        if cmd == 0:
            cnt += 3
            for k in range(cnt):
                outbuff[rpos + k] = inbuff[ipos]
            rpos += cnt
            ipos += 1

        # long RLE
        elif cmd == 1:
            cnt += <uint16_t>inbuff[ipos] << 4
            cnt += 19
            ipos += 1
            for k in range(cnt):
                outbuff[rpos + k] = inbuff[ipos]
            rpos += cnt
            ipos += 1

        # long pattern
        elif cmd == 2:
            ofs = cnt + 3
            ofs += <uint16_t>inbuff[ipos] << 4
            ipos += 1
            cnt = <uint16_t>inbuff[ipos]
            ipos += 1
            cnt += 16
            for k in range(cnt):
                outbuff[rpos + k] = outbuff[rpos - int(ofs) + k]
            rpos += cnt

        # short pattern
        elif (cmd >= 3) & (cmd <= 15):
            ofs = cnt + 3
            ofs += <uint16_t>inbuff[ipos] << 4
            ipos += 1
            for k in range(cmd):
                outbuff[rpos + k] = outbuff[rpos - int(ofs) + k]
            rpos += cmd

        else:
            raise ValueError("unknown RDC command")

    if len(outbuff) != result_length:
        raise ValueError("RDC: %v != %v\n", len(outbuff), result_length)

    return np.asarray(outbuff).tostring()

def process_byte_array_with_data(parser, int offset, int length, np.ndarray[uint8_t, ndim=2] byte_chunk,
                                 np.ndarray[dtype=object, ndim=2] string_chunk):

    cdef int s
    cdef int j
    cdef int m
    cdef int start
    cdef int end
    cdef bytes source
    cdef bytes temp
    cdef int jb
    cdef int js

    if (parser.compression != "") and (length < parser.row_length):
        source = parser._decompress(parser.row_length, parser._cached_page[offset:offset + length])
    else:
        source = parser._cached_page[offset:offset + length]

    s = 8 * parser._current_row_in_chunk_index
    js = 0
    jb = 0
    for j in range(parser.column_count):
        length = parser._column_data_lengths[j]
        if length == 0:
            break
        start = parser._column_data_offsets[j]
        end = start + length
        temp = source[start:end]
        if parser.column_types[j] == b'd':
            m = 8 - length
            if parser.byte_order == "<":
                byte_chunk[jb, s+m:s+8] = np.frombuffer(temp, dtype=np.uint8)
            else:
                byte_chunk[jb, s:s+length] = np.frombuffer(temp, dtype=np.uint8)
            jb += 1
        elif parser.column_types[j] == b's':
            string_chunk[js, parser._current_row_in_chunk_index] = bytes(temp)
            js += 1
        else:
            raise ValueError("unknown column type: %s" % parser.columns[j].ctype)

    parser._current_row_on_page_index += 1
    parser._current_row_in_chunk_index += 1
    parser._current_row_in_file_index += 1
