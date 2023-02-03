.. _copy_on_write::

{{ header }}

*************
Copy on write
*************

Copy on Write is a mechanism to simplify the indexing API and improve
performance through avoiding copies if possible.
CoW means that any DataFrame or Series derived from another in any way always
behaves as a copy.

Reference tracking
------------------

To be able to determine, if we have to make a copy when writing into a DataFrame,
we have to be aware, if the values are shared with another DataFrame. pandas
keeps track of all ``Blocks`` that share values with another block internally to
be able to tell when a copy needs to be triggered. The reference tracking
mechanism is implemented on the Block level.

We use a custom reference tracker object, e.g. ``BlockValuesRefs`` that keeps
track of every block, whose values share memory with each other. The reference
is held through a weak-reference. Every two blocks that share some memory should
point to the same ``BlockValuesRefs`` object. If one block goes out of
scope, the reference to this block dies. As a consequence, the reference tracker
object always knows how many blocks are alive and share memory.

We can ask the reference tracking object if there is another block alive that shares
data with us before writing into the values. We can trigger a copy before
writing if there is in fact another block alive.
