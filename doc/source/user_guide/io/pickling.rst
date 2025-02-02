.. _io.pickle:

========
Pickling
========

All pandas objects are equipped with ``to_pickle`` methods which use Python's
``cPickle`` module to save data structures to disk using the pickle format.

.. ipython:: python

   df
   df.to_pickle("foo.pkl")

The ``read_pickle`` function in the ``pandas`` namespace can be used to load
any pickled pandas object (or any other pickled object) from file:


.. ipython:: python

   pd.read_pickle("foo.pkl")

.. ipython:: python
   :suppress:

   os.remove("foo.pkl")

.. warning::

   Loading pickled data received from untrusted sources can be unsafe.

   See: https://docs.python.org/3/library/pickle.html

.. warning::

   :func:`read_pickle` is only guaranteed backwards compatible back to a few minor release.

.. _io.pickle.compression:

Compressed pickle files
-----------------------

:func:`read_pickle`, :meth:`DataFrame.to_pickle` and :meth:`Series.to_pickle` can read
and write compressed pickle files. The compression types of ``gzip``, ``bz2``, ``xz``, ``zstd`` are supported for reading and writing.
The ``zip`` file format only supports reading and must contain only one data file
to be read.

The compression type can be an explicit parameter or be inferred from the file extension.
If 'infer', then use ``gzip``, ``bz2``, ``zip``, ``xz``, ``zstd`` if filename ends in ``'.gz'``, ``'.bz2'``, ``'.zip'``,
``'.xz'``, or ``'.zst'``, respectively.

The compression parameter can also be a ``dict`` in order to pass options to the
compression protocol. It must have a ``'method'`` key set to the name
of the compression protocol, which must be one of
{``'zip'``, ``'gzip'``, ``'bz2'``, ``'xz'``, ``'zstd'``}. All other key-value pairs are passed to
the underlying compression library.

.. ipython:: python

   df = pd.DataFrame(
       {
           "A": np.random.randn(1000),
           "B": "foo",
           "C": pd.date_range("20130101", periods=1000, freq="s"),
       }
   )
   df

Using an explicit compression type:

.. ipython:: python

   df.to_pickle("data.pkl.compress", compression="gzip")
   rt = pd.read_pickle("data.pkl.compress", compression="gzip")
   rt

Inferring compression type from the extension:

.. ipython:: python

   df.to_pickle("data.pkl.xz", compression="infer")
   rt = pd.read_pickle("data.pkl.xz", compression="infer")
   rt

The default is to 'infer':

.. ipython:: python

   df.to_pickle("data.pkl.gz")
   rt = pd.read_pickle("data.pkl.gz")
   rt

   df["A"].to_pickle("s1.pkl.bz2")
   rt = pd.read_pickle("s1.pkl.bz2")
   rt

Passing options to the compression protocol in order to speed up compression:

.. ipython:: python

   df.to_pickle("data.pkl.gz", compression={"method": "gzip", "compresslevel": 1})

.. ipython:: python
   :suppress:

   os.remove("data.pkl.compress")
   os.remove("data.pkl.xz")
   os.remove("data.pkl.gz")
   os.remove("s1.pkl.bz2")

.. _io.msgpack:

msgpack
-------

pandas support for ``msgpack`` has been removed in version 1.0.0. It is
recommended to use :ref:`pickle <io.pickle>` instead.

Alternatively, you can also the Arrow IPC serialization format for on-the-wire
transmission of pandas objects. For documentation on pyarrow, see
`here <https://arrow.apache.org/docs/python/ipc.html>`__.
