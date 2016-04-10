import numpy as np
cimport numpy as np
from numpy cimport uint8_t, uint16_t, int8_t
import sas_constants as const

# rle_decompress decompresses data using a Run Length Encoding
# algorithm.  It is partially documented here:
#
# https://cran.r-project.org/web/packages/sas7bdat/vignettes/sas7bdat.pdf
cdef rle_decompress(int result_length, np.ndarray[uint8_t, ndim=1] inbuff):

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

    return np.asarray(result)


# rdc_decompress decompresses data using the Ross Data Compression algorithm:
#
#   http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/1992/9210/ross/ross.htm
cdef rdc_decompress(int result_length, np.ndarray[uint8_t, ndim=1] inbuff):

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

    return np.asarray(outbuff)

cdef decompress(object parser, int row_length, page):
    page = np.frombuffer(page, dtype=np.uint8)
    if parser.compression == const.rle_compression:
        return rle_decompress(row_length, page)
    elif parser.compression == const.rdc_compression:
        return rdc_decompress(row_length, page)
    else:
        raise ValueError("unknown SAS compression method: %s" %
                         parser.compression)


def do_read(object parser, int nrows):
    cdef int i

    for i in range(nrows):
        done = readline(parser)
        if done:
            break


cdef readline(object parser):

    cdef:
        int offset, bit_offset, align_correction, subheader_pointer_length

    bit_offset = parser._page_bit_offset
    subheader_pointer_length = parser._subheader_pointer_length

    # If there is no page, go to the end of the header and read a page.
    if parser._cached_page is None:
        parser._path_or_buf.seek(parser.header_length)
        done = parser._read_next_page()
        if done:
            return True

    # Loop until a data row is read
    while True:
        if parser._current_page_type == const.page_meta_type:
            flag = (parser._current_row_on_page_index >=
                    len(parser._current_page_data_subheader_pointers))
            if flag:
                done = parser._read_next_page()
                if done:
                    return True
                parser._current_row_on_page_index = 0
                continue
            current_subheader_pointer = (
                parser._current_page_data_subheader_pointers[
                    parser._current_row_on_page_index])
            process_byte_array_with_data(parser,
                                         current_subheader_pointer.offset,
                                         current_subheader_pointer.length)
            return False
        elif parser._current_page_type in const.page_mix_types:
            align_correction = (bit_offset + const.subheader_pointers_offset +
                                parser._current_page_subheaders_count *
                                subheader_pointer_length)
            align_correction = align_correction % 8
            offset = bit_offset + align_correction
            offset += const.subheader_pointers_offset
            offset += (parser._current_page_subheaders_count *
                       subheader_pointer_length)
            offset += parser._current_row_on_page_index * parser.row_length
            process_byte_array_with_data(parser, offset, parser.row_length)
            mn = min(parser.row_count, parser._mix_page_row_count)
            if parser._current_row_on_page_index == mn:
                done = parser._read_next_page()
                if done:
                    return True
                parser._current_row_on_page_index = 0
            return False
        elif parser._current_page_type == const.page_data_type:
            process_byte_array_with_data(parser,
                                         bit_offset +
                                         const.subheader_pointers_offset +
                                         parser._current_row_on_page_index *
                                         parser.row_length,
                                         parser.row_length)
            flag = (parser._current_row_on_page_index ==
                    parser._current_page_block_count)
            if flag:
                done = parser._read_next_page()
                if done:
                    return True
                parser._current_row_on_page_index = 0
            return False
        else:
            raise ValueError("unknown page type: %s",
                             parser._current_page_type)


cdef process_byte_array_with_data(object parser, int offset, int length):

    cdef:
        int s, j, k, m, start, jb, js, lngt
        long[:] lengths = parser._column_data_lengths
        long[:] offsets = parser._column_data_offsets
        char[:] column_types = parser.column_types
        uint8_t[:, :] byte_chunk = parser._byte_chunk
        object[:, :] string_chunk = parser._string_chunk
        np.ndarray[uint8_t, ndim=1] source
        np.ndarray[uint8_t, ndim=1] raw_source = np.frombuffer(parser._cached_page[offset:offset+length], dtype=np.uint8)

    if (parser.compression != "") and (length < parser.row_length):
        source = decompress(parser, parser.row_length, raw_source)
    else:
        source = raw_source

    s = 8 * parser._current_row_in_chunk_index
    js = 0
    jb = 0
    for j in range(parser.column_count):
        lngt = lengths[j]
        if lngt == 0:
            break
        start = offsets[j]
        if column_types[j] == b'd':
            if parser.byte_order == "<":
                m = s + 8 - lngt
            else:
                m = s
            for k in range(lngt):
                byte_chunk[jb, m + k] = source[start + k]
            jb += 1
        elif column_types[j] == b's':
            string_chunk[js, parser._current_row_in_chunk_index] = source[start:(start+lngt)].tostring().rstrip()
            js += 1
        else:
          raise ValueError("unknown column type: %s" % parser.columns[j].ctype)

    parser._current_row_on_page_index += 1
    parser._current_row_in_chunk_index += 1
    parser._current_row_in_file_index += 1
