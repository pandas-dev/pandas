.. _copy_on_write_dev:

{{ header }}

*************
Copy on write
*************

Copy on Write is a mechanism to simplify the indexing API and improve
performance through avoiding copies if possible.
CoW means that any DataFrame or Series derived from another in any way always
behaves as a copy. An explanation on how to use Copy on Write efficiently can be
found :ref:`here <copy_on_write>`.

Reference tracking
------------------

To be able to determine if we have to make a copy when writing into a DataFrame,
we have to be aware if the values are shared with another DataFrame. pandas
keeps track of all ``Blocks`` that share values with another block internally to
be able to tell when a copy needs to be triggered. The reference tracking
mechanism is implemented on the Block level.

We use a custom reference tracker object, ``BlockValuesRefs``, that keeps
track of every block, whose values share memory with each other. The reference
is held through a weak-reference. Every pair of blocks that share some memory should
point to the same ``BlockValuesRefs`` object. If one block goes out of
scope, the reference to this block dies. As a consequence, the reference tracker
object always knows how many blocks are alive and share memory.

Whenever a :class:`DataFrame` or :class:`Series` object is sharing data with another
object, it is required that each of those objects have its own BlockManager and Block
objects. Thus, in other words, one Block instance (that is held by a DataFrame, not
necessarily for intermediate objects) should always be uniquely used for only
a single DataFrame/Series object. For example, when you want to use the same
Block for another object, you can create a shallow copy of the Block instance
with ``block.copy(deep=False)`` (which will create a new Block instance with
the same underlying values and which will correctly set up the references).

We can ask the reference tracking object if there is another block alive that shares
data with us before writing into the values. We can trigger a copy before
writing if there is in fact another block alive.
