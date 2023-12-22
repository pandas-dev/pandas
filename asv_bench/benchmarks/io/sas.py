from pathlib import Path

from pandas import read_sas

ROOT = Path(__file__).parents[3] / "pandas" / "tests" / "io" / "sas" / "data"


class SAS:
    def time_read_sas7bdat(self):
        read_sas(ROOT / "test1.sas7bdat")

    def time_read_xpt(self):
        read_sas(ROOT / "paxraw_d_short.xpt")

    def time_read_sas7bdat_2(self):
        next(read_sas(ROOT / "0x00controlbyte.sas7bdat.bz2", chunksize=11000))

    def time_read_sas7bdat_2_chunked(self):
        for i, _ in enumerate(
            read_sas(ROOT / "0x00controlbyte.sas7bdat.bz2", chunksize=1000)
        ):
            if i == 10:
                break
