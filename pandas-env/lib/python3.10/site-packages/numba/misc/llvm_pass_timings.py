import re
import operator
import heapq
from collections import namedtuple
from collections.abc import Sequence
from contextlib import contextmanager
from functools import cached_property

from numba.core import config

import llvmlite.binding as llvm


class RecordLLVMPassTimings:
    """A helper context manager to track LLVM pass timings.
    """

    __slots__ = ["_data"]

    def __enter__(self):
        """Enables the pass timing in LLVM.
        """
        llvm.set_time_passes(True)
        return self

    def __exit__(self, exc_val, exc_type, exc_tb):
        """Reset timings and save report internally.
        """
        self._data = llvm.report_and_reset_timings()
        llvm.set_time_passes(False)
        return

    def get(self):
        """Retrieve timing data for processing.

        Returns
        -------
        timings: ProcessedPassTimings
        """
        return ProcessedPassTimings(self._data)


PassTimingRecord = namedtuple(
    "PassTimingRecord",
    [
        "user_time",
        "user_percent",
        "system_time",
        "system_percent",
        "user_system_time",
        "user_system_percent",
        "wall_time",
        "wall_percent",
        "pass_name",
        "instruction",
    ],
)


def _adjust_timings(records):
    """Adjust timing records because of truncated information.

    Details: The percent information can be used to improve the timing
    information.

    Returns
    -------
    res: List[PassTimingRecord]
    """
    total_rec = records[-1]
    assert total_rec.pass_name == "Total"  # guard for implementation error

    def make_adjuster(attr):
        time_attr = f"{attr}_time"
        percent_attr = f"{attr}_percent"
        time_getter = operator.attrgetter(time_attr)

        def adjust(d):
            """Compute percent x total_time = adjusted"""
            total = time_getter(total_rec)
            adjusted = total * d[percent_attr] * 0.01
            d[time_attr] = adjusted
            return d

        return adjust

    # Make adjustment functions for each field
    adj_fns = [
        make_adjuster(x) for x in ["user", "system", "user_system", "wall"]
    ]

    # Extract dictionaries from the namedtuples
    dicts = map(lambda x: x._asdict(), records)

    def chained(d):
        # Chain the adjustment functions
        for fn in adj_fns:
            d = fn(d)
        # Reconstruct the namedtuple
        return PassTimingRecord(**d)

    return list(map(chained, dicts))


class ProcessedPassTimings:
    """A class for processing raw timing report from LLVM.

    The processing is done lazily so we don't waste time processing unused
    timing information.
    """

    def __init__(self, raw_data):
        self._raw_data = raw_data

    def __bool__(self):
        return bool(self._raw_data)

    def get_raw_data(self):
        """Returns the raw string data.

        Returns
        -------
        res: str
        """
        return self._raw_data

    def get_total_time(self):
        """Compute the total time spend in all passes.

        Returns
        -------
        res: float
        """
        return self.list_records()[-1].wall_time

    def list_records(self):
        """Get the processed data for the timing report.

        Returns
        -------
        res: List[PassTimingRecord]
        """
        return self._processed

    def list_top(self, n):
        """Returns the top(n) most time-consuming (by wall-time) passes.

        Parameters
        ----------
        n: int
            This limits the maximum number of items to show.
            This function will show the ``n`` most time-consuming passes.

        Returns
        -------
        res: List[PassTimingRecord]
            Returns the top(n) most time-consuming passes in descending order.
        """
        records = self.list_records()
        key = operator.attrgetter("wall_time")
        return heapq.nlargest(n, records[:-1], key)

    def summary(self, topn=5, indent=0):
        """Return a string summarizing the timing information.

        Parameters
        ----------
        topn: int; optional
            This limits the maximum number of items to show.
            This function will show the ``topn`` most time-consuming passes.
        indent: int; optional
            Set the indentation level. Defaults to 0 for no indentation.

        Returns
        -------
        res: str
        """
        buf = []
        prefix = " " * indent

        def ap(arg):
            buf.append(f"{prefix}{arg}")

        ap(f"Total {self.get_total_time():.4f}s")
        ap("Top timings:")
        for p in self.list_top(topn):
            ap(f"  {p.wall_time:.4f}s ({p.wall_percent:5}%) {p.pass_name}")
        return "\n".join(buf)

    @cached_property
    def _processed(self):
        """A cached property for lazily processing the data and returning it.

        See ``_process()`` for details.
        """
        return self._process()

    def _process(self):
        """Parses the raw string data from LLVM timing report and attempts
        to improve the data by recomputing the times
        (See `_adjust_timings()``).
        """

        def parse(raw_data):
            """A generator that parses the raw_data line-by-line to extract
            timing information for each pass.
            """
            lines = raw_data.splitlines()
            colheader = r"[a-zA-Z+ ]+"
            # Take at least one column header.
            multicolheaders = fr"(?:\s*-+{colheader}-+)+"

            line_iter = iter(lines)
            # find column headers
            header_map = {
                "User Time": "user",
                "System Time": "system",
                "User+System": "user_system",
                "Wall Time": "wall",
                "Instr": "instruction",
                "Name": "pass_name",
            }
            for ln in line_iter:
                m = re.match(multicolheaders, ln)
                if m:
                    # Get all the column headers
                    raw_headers = re.findall(r"[a-zA-Z][a-zA-Z+ ]+", ln)
                    headers = [header_map[k.strip()] for k in raw_headers]
                    break

            assert headers[-1] == 'pass_name'
            # compute the list of available attributes from the column headers
            attrs = []
            n = r"\s*((?:[0-9]+\.)?[0-9]+)"
            pat = ""
            for k in headers[:-1]:
                if k == "instruction":
                    pat += n
                else:
                    attrs.append(f"{k}_time")
                    attrs.append(f"{k}_percent")
                    pat += rf"\s+(?:{n}\s*\({n}%\)|-+)"

            # put default value 0.0 to all missing attributes
            missing = {}
            for k in PassTimingRecord._fields:
                if k not in attrs and k != 'pass_name':
                    missing[k] = 0.0
            # parse timings
            pat += r"\s*(.*)"
            for ln in line_iter:
                m = re.match(pat, ln)
                if m is not None:
                    raw_data = list(m.groups())
                    data = {k: float(v) if v is not None else 0.0
                            for k, v in zip(attrs, raw_data)}
                    data.update(missing)
                    pass_name = raw_data[-1]
                    rec = PassTimingRecord(
                        pass_name=pass_name, **data,
                    )
                    yield rec
                    if rec.pass_name == "Total":
                        # "Total" means the report has ended
                        break
            # Check that we have reach the end of the report
            remaining = '\n'.join(line_iter)
            if remaining:
                raise ValueError(
                    f"unexpected text after parser finished:\n{remaining}"
                )

        # Parse raw data
        records = list(parse(self._raw_data))
        return _adjust_timings(records)


NamedTimings = namedtuple("NamedTimings", ["name", "timings"])


class PassTimingsCollection(Sequence):
    """A collection of pass timings.

    This class implements the ``Sequence`` protocol for accessing the
    individual timing records.
    """

    def __init__(self, name):
        self._name = name
        self._records = []

    @contextmanager
    def record(self, name):
        """Record new timings and append to this collection.

        Note: this is mainly for internal use inside the compiler pipeline.

        See also ``RecordLLVMPassTimings``

        Parameters
        ----------
        name: str
            Name for the records.
        """
        if config.LLVM_PASS_TIMINGS:
            # Recording of pass timings is enabled
            with RecordLLVMPassTimings() as timings:
                yield
            rec = timings.get()
            # Only keep non-empty records
            if rec:
                self._append(name, rec)
        else:
            # Do nothing. Recording of pass timings is disabled.
            yield

    def _append(self, name, timings):
        """Append timing records

        Parameters
        ----------
        name: str
            Name for the records.
        timings: ProcessedPassTimings
            the timing records.
        """
        self._records.append(NamedTimings(name, timings))

    def get_total_time(self):
        """Computes the sum of the total time across all contained timings.

        Returns
        -------
        res: float or None
            Returns the total number of seconds or None if no timings were
            recorded
        """
        if self._records:
            return sum(r.timings.get_total_time() for r in self._records)
        else:
            return None

    def list_longest_first(self):
        """Returns the timings in descending order of total time duration.

        Returns
        -------
        res: List[ProcessedPassTimings]
        """
        return sorted(self._records,
                      key=lambda x: x.timings.get_total_time(),
                      reverse=True)

    @property
    def is_empty(self):
        """
        """
        return not self._records

    def summary(self, topn=5):
        """Return a string representing the summary of the timings.

        Parameters
        ----------
        topn: int; optional, default=5.
            This limits the maximum number of items to show.
            This function will show the ``topn`` most time-consuming passes.

        Returns
        -------
        res: str

        See also ``ProcessedPassTimings.summary()``
        """
        if self.is_empty:
            return "No pass timings were recorded"
        else:
            buf = []
            ap = buf.append
            ap(f"Printing pass timings for {self._name}")
            overall_time = self.get_total_time()
            ap(f"Total time: {overall_time:.4f}")
            for i, r in enumerate(self._records):
                ap(f"== #{i} {r.name}")
                percent = r.timings.get_total_time() / overall_time * 100
                ap(f" Percent: {percent:.1f}%")
                ap(r.timings.summary(topn=topn, indent=1))
            return "\n".join(buf)

    def __getitem__(self, i):
        """Get the i-th timing record.

        Returns
        -------
        res: (name, timings)
            A named tuple with two fields:

            - name: str
            - timings: ProcessedPassTimings
        """
        return self._records[i]

    def __len__(self):
        """Length of this collection.
        """
        return len(self._records)

    def __str__(self):
        return self.summary()
