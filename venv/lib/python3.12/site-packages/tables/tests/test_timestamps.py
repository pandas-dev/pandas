"""This test unit checks control of dataset timestamps with track_times.

"""

import sys
import time
import hashlib
from pathlib import Path

import tables as tb
from tables.tests import common

HEXDIGEST = "2aafb84ab739bb4ae61d2939dc010bfd"


class Record(tb.IsDescription):
    var1 = tb.StringCol(itemsize=4)  # 4-character String
    var2 = tb.IntCol()  # integer
    var3 = tb.Int16Col()  # short integer


class TrackTimesMixin:
    def _add_datasets(self, group, j, track_times):
        # Create a table
        table = self.h5file.create_table(
            group,
            f"table{j}",
            Record,
            title=self.title,
            filters=None,
            track_times=track_times,
        )
        # Get the record object associated with the new table
        d = table.row
        # Fill the table
        for i in range(self.nrows):
            d["var1"] = "%04d" % (self.nrows - i)
            d["var2"] = i
            d["var3"] = i * 2
            d.append()  # This injects the Record values
        # Flush the buffer for this table
        table.flush()

        # Create a couple of arrays in each group
        var1List = [x["var1"] for x in table.iterrows()]
        var3List = [x["var3"] for x in table.iterrows()]

        self.h5file.create_array(
            group, f"array{j}", var1List, f"col {j}", track_times=track_times
        )

        # Create CArrays as well
        self.h5file.create_carray(
            group,
            name=f"carray{j}",
            obj=var3List,
            title=f"col {j + 2}",
            track_times=track_times,
        )

        # Create EArrays as well
        ea = self.h5file.create_earray(
            group,
            f"earray{j}",
            tb.StringAtom(itemsize=4),
            (0,),
            f"col {j + 4}",
            track_times=track_times,
        )
        # And fill them with some values
        ea.append(var1List)

        # Finally VLArrays too
        vla = self.h5file.create_vlarray(
            group,
            f"vlarray{j}",
            tb.Int16Atom(),
            f"col {j + 6}",
            track_times=track_times,
        )
        # And fill them with some values
        vla.append(var3List)


class TimestampTestCase(
    TrackTimesMixin, common.TempFileMixin, common.PyTablesTestCase
):
    title = "A title"
    nrows = 10

    def setUp(self):
        super().setUp()
        self.populateFile()

    def populateFile(self):
        group = self.h5file.root
        for j in range(4):
            track_times = bool(j % 2)
            self._add_datasets(group, j, track_times)

    def test00_checkTimestamps(self):
        """Checking retrieval of timestamps"""

        for pattern in (
            "/table{}",
            "/array{}",
            "/carray{}",
            "/earray{}",
            "/vlarray{}",
        ):
            # Verify that:
            # - if track_times was False, ctime is 0
            # - if track_times was True, ctime is not 0,
            #   and has either stayed the same or incremented
            tracked_ctimes = []
            for j in range(4):
                track_times = bool(j % 2)
                node = pattern.format(j)
                obj = self.h5file.get_node(node)
                # Test property retrieval
                self.assertEqual(obj.track_times, track_times)
                timestamps = obj._get_obj_timestamps()
                self.assertEqual(timestamps.atime, 0)
                self.assertEqual(timestamps.mtime, 0)
                self.assertEqual(timestamps.btime, 0)
                if not track_times:
                    self.assertEqual(timestamps.ctime, 0)
                else:
                    self.assertNotEqual(timestamps.ctime, 0)
                    tracked_ctimes.append(timestamps.ctime)
            self.assertGreaterEqual(tracked_ctimes[1], tracked_ctimes[0])


class BitForBitTestCase(
    TrackTimesMixin, common.TempFileMixin, common.PyTablesTestCase
):
    title = "A title"
    nrows = 10

    def repopulateFile(self, track_times):
        self.h5file.close()
        self.h5file = tb.open_file(self.h5fname, mode="w")
        group = self.h5file.root
        self._add_datasets(group, 1, track_times)
        self.h5file.close()

    def test00_checkReproducibility(self):
        """Checking bit-for-bit reproducibility with no track_times"""

        self.repopulateFile(track_times=False)
        hexdigest_wo_track_1 = self._get_digest(self.h5fname)
        self.repopulateFile(track_times=True)
        hexdigest_w_track_1 = self._get_digest(self.h5fname)
        time.sleep(1)
        self.repopulateFile(track_times=True)
        hexdigest_w_track_2 = self._get_digest(self.h5fname)
        self.repopulateFile(track_times=False)
        hexdigest_wo_track_2 = self._get_digest(self.h5fname)
        self.assertEqual(HEXDIGEST, hexdigest_wo_track_1)
        self.assertEqual(hexdigest_wo_track_1, hexdigest_wo_track_2)
        self.assertNotEqual(hexdigest_wo_track_1, hexdigest_w_track_1)
        self.assertNotEqual(hexdigest_w_track_1, hexdigest_w_track_2)

    def _get_digest(self, filename):
        md5 = hashlib.md5()
        for data in Path(filename).read_bytes():
            md5.update(data)

        hexdigest = md5.hexdigest()

        return hexdigest


def suite():
    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1 # Uncomment this only for testing purposes!

    for i in range(niter):
        theSuite.addTest(common.make_suite(TimestampTestCase))
        theSuite.addTest(common.make_suite(BitForBitTestCase))

    return theSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
