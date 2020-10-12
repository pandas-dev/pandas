.. _timedeltas:

{{ header }}

.. _timedeltas.timedeltas:

***********
Time deltas
***********

Timedeltas are differences in times, expressed in difference units, e.g. days, hours, minutes,
seconds. They can be both positive and negative.

``Timedelta`` is a subclass of ``datetime.timedelta``, and behaves in a similar manner,
but allows compatibility with ``np.timedelta64`` types as well as a host of custom representation,
parsing, and attributes.

Parsing
-------

You can construct a ``Timedelta`` scalar through various arguments, including `ISO 8601 Duration`_ strings.

.. ipython:: python

   import datetime

   # strings
   pd.Timedelta("1 days")
   pd.Timedelta("1 days 00:00:00")
   pd.Timedelta("1 days 2 hours")
   pd.Timedelta("-1 days 2 min 3us")

   # like datetime.timedelta
   # note: these MUST be specified as keyword arguments
   pd.Timedelta(days=1, seconds=1)

   # integers with a unit
   pd.Timedelta(1, unit="d")

   # from a datetime.timedelta/np.timedelta64
   pd.Timedelta(datetime.timedelta(days=1, seconds=1))
   pd.Timedelta(np.timedelta64(1, "ms"))

   # negative Timedeltas have this string repr
   # to be more consistent with datetime.timedelta conventions
   pd.Timedelta("-1us")

   # a NaT
   pd.Timedelta("nan")
   pd.Timedelta("nat")

   # ISO 8601 Duration strings
   pd.Timedelta("P0DT0H1M0S")
   pd.Timedelta("P0DT0H0M0.000000123S")

:ref:`DateOffsets<timeseries.offsets>` (``Day, Hour, Minute, Second, Milli, Micro, Nano``) can also be used in construction.

.. ipython:: python

   pd.Timedelta(pd.offsets.Second(2))

Further, operations among the scalars yield another scalar ``Timedelta``.

.. ipython:: python

   pd.Timedelta(pd.offsets.Day(2)) + pd.Timedelta(pd.offsets.Second(2)) + pd.Timedelta(
       "00:00:00.000123"
   )

to_timedelta
~~~~~~~~~~~~

Using the top-level ``pd.to_timedelta``, you can convert a scalar, array, list,
or Series from a recognized timedelta format / value into a ``Timedelta`` type.
It will construct Series if the input is a Series, a scalar if the input is
scalar-like, otherwise it will output a ``TimedeltaIndex``.

You can parse a single string to a Timedelta:

.. ipython:: python

   pd.to_timedelta("1 days 06:05:01.00003")
   pd.to_timedelta("15.5us")

or a list/array of strings:

.. ipython:: python

   pd.to_timedelta(["1 days 06:05:01.00003", "15.5us", "nan"])

The ``unit`` keyword argument specifies the unit of the Timedelta:

.. ipython:: python

   pd.to_timedelta(np.arange(5), unit="s")
   pd.to_timedelta(np.arange(5), unit="d")

.. _timedeltas.limitations:

Timedelta limitations
~~~~~~~~~~~~~~~~~~~~~

pandas represents ``Timedeltas`` in nanosecond resolution using
64 bit integers. As such, the 64 bit integer limits determine
the ``Timedelta`` limits.

.. ipython:: python

   pd.Timedelta.min
   pd.Timedelta.max

.. _timedeltas.operations:

Operations
----------

You can operate on Series/DataFrames and construct ``timedelta64[ns]`` Series through
subtraction operations on ``datetime64[ns]`` Series, or ``Timestamps``.

.. ipython:: python

   s = pd.Series(pd.date_range("2012-1-1", periods=3, freq="D"))
   td = pd.Series([pd.Timedelta(days=i) for i in range(3)])
   df = pd.DataFrame({"A": s, "B": td})
   df
   df["C"] = df["A"] + df["B"]
   df
   df.dtypes

   s - s.max()
   s - datetime.datetime(2011, 1, 1, 3, 5)
   s + datetime.timedelta(minutes=5)
   s + pd.offsets.Minute(5)
   s + pd.offsets.Minute(5) + pd.offsets.Milli(5)

Operations with scalars from a ``timedelta64[ns]`` series:

.. ipython:: python

   y = s - s[0]
   y

Series of timedeltas with ``NaT`` values are supported:

.. ipython:: python

   y = s - s.shift()
   y

Elements can be set to ``NaT`` using ``np.nan`` analogously to datetimes:

.. ipython:: python

   y[1] = np.nan
   y

Operands can also appear in a reversed order (a singular object operated with a Series):

.. ipython:: python

   s.max() - s
   datetime.datetime(2011, 1, 1, 3, 5) - s
   datetime.timedelta(minutes=5) + s

``min, max`` and the corresponding ``idxmin, idxmax`` operations are supported on frames:

.. ipython:: python

   A = s - pd.Timestamp("20120101") - pd.Timedelta("00:05:05")
   B = s - pd.Series(pd.date_range("2012-1-2", periods=3, freq="D"))

   df = pd.DataFrame({"A": A, "B": B})
   df

   df.min()
   df.min(axis=1)

   df.idxmin()
   df.idxmax()

``min, max, idxmin, idxmax`` operations are supported on Series as well. A scalar result will be a ``Timedelta``.

.. ipython:: python

   df.min().max()
   df.min(axis=1).min()

   df.min().idxmax()
   df.min(axis=1).idxmin()

You can fillna on timedeltas, passing a timedelta to get a particular value.

.. ipython:: python

   y.fillna(pd.Timedelta(0))
   y.fillna(pd.Timedelta(10, unit="s"))
   y.fillna(pd.Timedelta("-1 days, 00:00:05"))

You can also negate, multiply and use ``abs`` on ``Timedeltas``:

.. ipython:: python

   td1 = pd.Timedelta("-1 days 2 hours 3 seconds")
   td1
   -1 * td1
   -td1
   abs(td1)

.. _timedeltas.timedeltas_reductions:

Reductions
----------

Numeric reduction operation for ``timedelta64[ns]`` will return ``Timedelta`` objects. As usual
``NaT`` are skipped during evaluation.

.. ipython:: python

   y2 = pd.Series(
       pd.to_timedelta(["-1 days +00:00:05", "nat", "-1 days +00:00:05", "1 days"])
   )
   y2
   y2.mean()
   y2.median()
   y2.quantile(0.1)
   y2.sum()

.. _timedeltas.timedeltas_convert:

Frequency conversion
--------------------

Timedelta Series, ``TimedeltaIndex``, and ``Timedelta`` scalars can be converted to other 'frequencies' by dividing by another timedelta,
or by astyping to a specific timedelta type. These operations yield Series and propagate ``NaT`` -> ``nan``.
Note that division by the NumPy scalar is true division, while astyping is equivalent of floor division.

.. ipython:: python

   december = pd.Series(pd.date_range("20121201", periods=4))
   january = pd.Series(pd.date_range("20130101", periods=4))
   td = january - december

   td[2] += datetime.timedelta(minutes=5, seconds=3)
   td[3] = np.nan
   td

   # to days
   td / np.timedelta64(1, "D")
   td.astype("timedelta64[D]")

   # to seconds
   td / np.timedelta64(1, "s")
   td.astype("timedelta64[s]")

   # to months (these are constant months)
   td / np.timedelta64(1, "M")

Dividing or multiplying a ``timedelta64[ns]`` Series by an integer or integer Series
yields another ``timedelta64[ns]`` dtypes Series.

.. ipython:: python

   td * -1
   td * pd.Series([1, 2, 3, 4])

Rounded division (floor-division) of a ``timedelta64[ns]`` Series by a scalar
``Timedelta`` gives a series of integers.

.. ipython:: python

   td // pd.Timedelta(days=3, hours=4)
   pd.Timedelta(days=3, hours=4) // td

.. _timedeltas.mod_divmod:

The mod (%) and divmod operations are defined for ``Timedelta`` when operating with another timedelta-like or with a numeric argument.

.. ipython:: python

   pd.Timedelta(hours=37) % datetime.timedelta(hours=2)

   # divmod against a timedelta-like returns a pair (int, Timedelta)
   divmod(datetime.timedelta(hours=2), pd.Timedelta(minutes=11))

   # divmod against a numeric returns a pair (Timedelta, Timedelta)
   divmod(pd.Timedelta(hours=25), 86400000000000)

Attributes
----------

You can access various components of the ``Timedelta`` or ``TimedeltaIndex`` directly using the attributes ``days,seconds,microseconds,nanoseconds``. These are identical to the values returned by ``datetime.timedelta``, in that, for example, the ``.seconds`` attribute represents the number of seconds >= 0 and < 1 day. These are signed according to whether the ``Timedelta`` is signed.

These operations can also be directly accessed via the ``.dt`` property of the ``Series`` as well.

.. note::

   Note that the attributes are NOT the displayed values of the ``Timedelta``. Use ``.components`` to retrieve the displayed values.

For a ``Series``:

.. ipython:: python

   td.dt.days
   td.dt.seconds

You can access the value of the fields for a scalar ``Timedelta`` directly.

.. ipython:: python

   tds = pd.Timedelta("31 days 5 min 3 sec")
   tds.days
   tds.seconds
   (-tds).seconds

You can use the ``.components`` property to access a reduced form of the timedelta. This returns a ``DataFrame`` indexed
similarly to the ``Series``. These are the *displayed* values of the ``Timedelta``.

.. ipython:: python

   td.dt.components
   td.dt.components.seconds

.. _timedeltas.isoformat:

You can convert a ``Timedelta`` to an `ISO 8601 Duration`_ string with the
``.isoformat`` method

.. ipython:: python

    pd.Timedelta(
        days=6, minutes=50, seconds=3, milliseconds=10, microseconds=10, nanoseconds=12
    ).isoformat()

.. _ISO 8601 Duration: https://en.wikipedia.org/wiki/ISO_8601#Durations

.. _timedeltas.index:

TimedeltaIndex
--------------

To generate an index with time delta, you can use either the :class:`TimedeltaIndex` or
the :func:`timedelta_range` constructor.

Using ``TimedeltaIndex`` you can pass string-like, ``Timedelta``, ``timedelta``,
or ``np.timedelta64`` objects. Passing ``np.nan/pd.NaT/nat`` will represent missing values.

.. ipython:: python

   pd.TimedeltaIndex(
       [
           "1 days",
           "1 days, 00:00:05",
           np.timedelta64(2, "D"),
           datetime.timedelta(days=2, seconds=2),
       ]
   )

The string 'infer' can be passed in order to set the frequency of the index as the
inferred frequency upon creation:

.. ipython:: python

   pd.TimedeltaIndex(["0 days", "10 days", "20 days"], freq="infer")

Generating ranges of time deltas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to :func:`date_range`, you can construct regular ranges of a ``TimedeltaIndex``
using :func:`timedelta_range`.  The default frequency for ``timedelta_range`` is
calendar day:

.. ipython:: python

   pd.timedelta_range(start="1 days", periods=5)

Various combinations of ``start``, ``end``, and ``periods`` can be used with
``timedelta_range``:

.. ipython:: python

   pd.timedelta_range(start="1 days", end="5 days")

   pd.timedelta_range(end="10 days", periods=4)

The ``freq`` parameter can passed a variety of :ref:`frequency aliases <timeseries.offset_aliases>`:

.. ipython:: python

   pd.timedelta_range(start="1 days", end="2 days", freq="30T")

   pd.timedelta_range(start="1 days", periods=5, freq="2D5H")


Specifying ``start``, ``end``, and ``periods`` will generate a range of evenly spaced
timedeltas from ``start`` to ``end`` inclusively, with ``periods`` number of elements
in the resulting ``TimedeltaIndex``:

.. ipython:: python

   pd.timedelta_range("0 days", "4 days", periods=5)

   pd.timedelta_range("0 days", "4 days", periods=10)

Using the TimedeltaIndex
~~~~~~~~~~~~~~~~~~~~~~~~

Similarly to other of the datetime-like indices, ``DatetimeIndex`` and ``PeriodIndex``, you can use
``TimedeltaIndex`` as the index of pandas objects.

.. ipython:: python

   s = pd.Series(
       np.arange(100),
       index=pd.timedelta_range("1 days", periods=100, freq="h"),
   )
   s

Selections work similarly, with coercion on string-likes and slices:

.. ipython:: python

   s["1 day":"2 day"]
   s["1 day 01:00:00"]
   s[pd.Timedelta("1 day 1h")]

Furthermore you can use partial string selection and the range will be inferred:

.. ipython:: python

   s["1 day":"1 day 5 hours"]

Operations
~~~~~~~~~~

Finally, the combination of ``TimedeltaIndex`` with ``DatetimeIndex`` allow certain combination operations that are NaT preserving:

.. ipython:: python

   tdi = pd.TimedeltaIndex(["1 days", pd.NaT, "2 days"])
   tdi.to_list()
   dti = pd.date_range("20130101", periods=3)
   dti.to_list()
   (dti + tdi).to_list()
   (dti - tdi).to_list()

Conversions
~~~~~~~~~~~

Similarly to frequency conversion on a ``Series`` above, you can convert these indices to yield another Index.

.. ipython:: python

   tdi / np.timedelta64(1, "s")
   tdi.astype("timedelta64[s]")

Scalars type ops work as well. These can potentially return a *different* type of index.

.. ipython:: python

   # adding or timedelta and date -> datelike
   tdi + pd.Timestamp("20130101")

   # subtraction of a date and a timedelta -> datelike
   # note that trying to subtract a date from a Timedelta will raise an exception
   (pd.Timestamp("20130101") - tdi).to_list()

   # timedelta + timedelta -> timedelta
   tdi + pd.Timedelta("10 days")

   # division can result in a Timedelta if the divisor is an integer
   tdi / 2

   # or a Float64Index if the divisor is a Timedelta
   tdi / tdi[0]

.. _timedeltas.resampling:

Resampling
----------

Similar to :ref:`timeseries resampling <timeseries.resampling>`, we can resample with a ``TimedeltaIndex``.

.. ipython:: python

   s.resample("D").mean()
