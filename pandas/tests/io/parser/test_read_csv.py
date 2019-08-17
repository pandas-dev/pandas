import pytest
import csv
import os
import sys


rootpath = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(rootpath, "pandas/io"))
from pandas.io import parsers



TMP_PATH = "tmp"

"""
To run test, run 'python path/to/test_parsers.py'
"""


"""
 test_read_csv_without_encoding_kwarg returns result of read_csv method.
    - if exception is raised from method, the result returned is the exception.
"""
def test_read_csv_without_encoding_kwarg(file):
    try:
        result = parsers.read_csv(file)
    except Exception as e:
        result = e
    return result


def write_csv_file(filename, data, encoding, delimiter=",", newline=""):
    with open(filename, 'w', newline=newline, encoding=encoding) as csv_file:
        writer = csv.writer(csv_file, delimiter=delimiter)
        # for row in data:
        writer.writerow(data)
    return filename


def test():
    test_results = {}
    test_dtypes = ['ascii',
                     'big5',
                     'big5hkscs',
                     'cp037',
                     'cp273',
                     'cp424',
                     'cp437',
                     'cp500',
                     'cp720',
                     'cp737',
                     'cp775',
                     'cp850',
                     'cp852',
                     'cp855',
                     'cp856',
                     'cp857',
                     'cp858',
                     'cp860',
                     'cp861',
                     'cp862',
                     'cp863',
                     'cp864',
                     'cp865',
                     'cp866',
                     'cp869',
                     'cp874',
                     'cp875',
                     'cp932',
                     'cp949',
                     'cp950',
                     'cp1006',
                     'cp1026',
                     'cp1125',
                     'cp1140',
                     'cp1250',
                     'cp1251',
                     'cp1252',
                     'cp1253',
                     'cp1254',
                     'cp1255',
                     'cp1256',
                     'cp1257',
                     'cp1258',
                     'cp65001',
                     'euc_jp',
                     'euc_jis_2004',
                     'euc_jisx0213',
                     'euc_kr',
                     'gb2312',
                     'gbk',
                     'gb18030',
                     'hz',
                     'iso2022_jp',
                     'iso2022_jp_1',
                     'iso2022_jp_2',
                     'iso2022_jp_2004',
                     'iso2022_jp_3',
                     'iso2022_jp_ext',
                     'iso2022_kr',
                     'latin_1',
                     'iso8859_2',
                     'iso8859_3',
                     'iso8859_4',
                     'iso8859_5',
                     'iso8859_6',
                     'iso8859_7',
                     'iso8859_8',
                     'iso8859_9',
                     'iso8859_10',
                     'iso8859_11',
                     'iso8859_13',
                     'iso8859_14',
                     'iso8859_15',
                     'iso8859_16',
                     'johab',
                     'koi8_r',
                     'koi8_t',
                     'koi8_u',
                     'kz1048',
                     'mac_cyrillic',
                     'mac_greek',
                     'mac_iceland',
                     'mac_latin2',
                     'mac_roman',
                     'mac_turkish',
                     'ptcp154',
                     'shift_jis',
                     'shift_jis_2004',
                     'shift_jisx0213',
                     'utf_32',
                     'utf_32_be',
                     'utf_32_le',
                     'utf_16',
                     'utf_16_be',
                     'utf_16_le',
                     'utf_7',
                     'utf_8',
                     'utf_8_sig'
                     ]
    data = """
            one,two,three
            1,2,3
            uno,dos,tres
            """
    for i, dtype in enumerate(test_dtypes):
        file = write_csv_file(f"test{i}.csv", data, dtype)
        result = test_read_csv_without_encoding_kwarg(file)
        test_results[dtype] = result

    print("test results: ", test_results)


if __name__ == '__main__':
    test()
