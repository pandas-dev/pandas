.. _roadmap.indexing_views:

==================
Indexing and Views
==================

*A proposal for consistent, clear copy vs. view semantics in pandas' indexing.*

**Issue**: https://github.com/pandas-dev/pandas/issues/36195

Motivation
----------

pandas’ current behavior on whether indexing returns a view or copy is
confusing. Even for experienced users, it’s hard to tell whether a view or copy
will be returned (see below for a summary). We’d like to provide an API that is
consistent and sensible about returning views vs. copies.

We also care about performance. Returning views from indexing operations is
faster and reduces memory usage (at least for that operation; whether it’s
faster for a full workflow depends on whether downstream operations trigger a
copy (possibly through block consolidation)).

Finally, there are API / usability issues around views. It can be challenging to
know the user’s intent in operations that modify a subset of a DataFrame (column
and/or row selection), like:

.. code-block:: python

   >>> df = pd.DataFrame({"A”": [1, 2], "B": [3, 4]})
   >>> df2 = df[["A"]]
   >>> df2.iloc[:, 0] = 10

Did the user intend to modify ``df`` when they modified ``df2`` (setting aside
issues with the current implementation)? In other words, if we had a perfectly
consistent world where indexing the columns always returned views or always
returned a copy, does the code above imply that the user wants to mutate ``df``?

There are two possible behaviours the user might intend:

1. I know my subset might be a view of the original and I want to modify the
   original as well.
2. I just want to modify the subset without modifying the original.

Today, pandas’ inconsistency means neither of these workflows is really
possible. The first is difficult, because indexing operations often (though not
always) return copies, and even when a view is returned you sometimes get a
``SettingWithCopyWarning`` when mutating. The second is somewhat possible, but
requires many defensive copies (to avoid ``SettingWithCopyWarning``, or to
ensure that you have a copy when a view was returned).

Proposal Summary
----------------

For these reasons (consistency, performance, code clarity), we propose three
changes:

1. Indexing always returns a view when possible. This means that indexing
   columns of a dataframe always returns a view
   (https://github.com/pandas-dev/pandas/pull/33597), and indexing rows may
   return a view, depending on the type of the row indexer.
2. We implement Error-on-Write (explained below)
3. We provide APIs for explicitly marking a DataFrame as a “mutable view”
   (mutating the dataframe would mutate its parents) and copying a dataframe
   only if needed to avoid concerns with mutating other dataframes (i.e. it is
   not a view on another dataframe).

The intent is to capture the performance benefits of views, while allowing users
to explicitly choose the behavior they want for inplace operations that might
mutate other dataframes. This essentially makes returning views an internal
optimization, without the user needing to know if the specific indexing
operation would return a view or a copy.

Taking the example from above, if the user wants to make use of the fact that
``df2`` is a view to modify the original ``df``, they would write:

.. code-block:: python

   # Case 1: user wants mutations of df2 to be reflected in df
   >>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
   >>> df2 = df[["A"]].as_mutable_view() # name TBD
   >>> df2.iloc[:, 0] = 10
   >>> df.iloc[0, 0] # df was mutated 10

For the user who wishes to not mutate the parent, we require that the user
explicitly break the reference from ``df2`` to ``df`` by implementing “Error on Write”.

.. code-block:: python

   # Case 2: The user does not want mutating df2 to mutate df, via EoW
   >>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
   >>> df2 = df[["A"]]
   >>> df2.iloc[0, 0] = 10
   MutableViewError("error on write to subset of other dataframe")
   >>> df2 = df2.copy_if_needed() # API is TBD. Could be a keyword argument to copy.
   >>> df2.iloc[:, 0] = 10
   >>> df.iloc[0, 0] # df was not mutated 1

Copy-on-Write vs. Error-on-Write
--------------------------------

Consider the following example:

.. code-block:: python

   >>> df2 = df[['A']]
   >>> df2.iloc[0, 0] = 10 # df2 can be a view of df, what happens by default?
   >>> df3 = df[df['A'] == 1]
   >>> df3.iloc[0, 0] = 10 # df3 is already a copy of df, what happens by default?

We have a few options for the default:

1. Well-Defined copy/view rules: ensure we have more consistent rules (e.g.
   selecting columns is always a view), and then views result in mutating the
   parent, copies not. This comes down to fixing some bugs and clearly
   documenting and testing which operations are views, and which are copies.
2. Copy-on-Write: The setitem would check if it’s a view on another dataframe.
   If it is, then we would copy our data before mutating.
3. Error-on-Write: The setitem would check if it’s a subset of another dataframe
   (both view of copy). Only rather than copying in case of a view we would
   raise an exception telling the user to either copy the data with
   ``.copy_if_needed()`` (name TBD) or mark the frame as “a mutable view” with
   ``.as_mutable_view()`` (name TBD).

We propose "Error on Write" by default. This forces a decision on the user, and
is the most explicit in terms of code.

Additionally, consider the "classic" case of chained indexing, which was the
original motivation for the ``SettingWithCopy`` warning

.. code-block:: python

   >>> df[df['B'] > 4]['B'] = 10

That is roughly equivalent to

.. code-block:: python

   >>> df2 = df[df['B'] > 4] # Copy under NumPy’s rules
   >>> df2['B'] = 10 # Update (the copy) df2, df not changed
   >>> del df2 # All references to df2 are lost, goes out of scope

And so ``df`` is not modified. If we adopted Copy On Write to completely replace the
current ``SettingWithCopy`` warning, we would restore the old behavior of silently
“failing” to update ``df2``. Under Error on Write, we’d track that the ``df2`` created
by the first getitem references ``df`` and raise an exception when it was being
mutated.

New methods
-----------

In addition to the behavior changes to indexing columns, this proposal includes
two new methods for controlling behavior in operations downstream of an indexing
operation.

.. code-block:: python

   def as_mutable_view(self):  # name TBD
       """
       Mark a DataFrame as mutable so that setitem operations propagate.
   
       Any setitem operations on the returned DataFrame will propagate
       to the DataFrame(s) this DataFrame is a view on.
   
       Examples
       --------
       >>> df1 = pd.DataFrame({"A": [1, 2]})
       >>> df2 = df[["A"]].as_mutable_view()  # df2 is a view on df
       >>> df2.iloc[0, 0] = 10
       >>> df1.iloc[0, 0]  # The parent df1 was mutated.
       10
       """

If we implement Error-On-Write, a ``copy_if_needed`` method is necessary for
libraries and user code to avoid unnecessary defensive copying.

.. code-block:: python

   def copy_if_needed(self):  # name TBD
       """
       Copy the data in a Series / DataFrame if it is a view on some other.
   
       This will copy the data backing a DataFrame only if it's a view
       on other some other dataframe. If it's not a view then no data is
       copied.
   
       Examples
       --------
       >>> df1 = pd.DataFrame({"A": [1, 2]})
       >>> df2 = df1[["A"]]  # df2 is a view on df1
       >>> df3 = df2.copy_if_needed()  # triggers a copy
   
       When no copy is necessary (the object is not a view on another dataframe)
       then no copy is performed.
   
       >>> df4 = df1[df1['a'] == 1].copy_if_needed()  # No copy, since boolean masking already returned a copy
       """


These two methods give users the control to say whether setitem operations on a
dataframe that is a view on another dataframe should mutate the “parent”
dataframe. Users wishing to mutate the parent will make it explicit with
``.as_mutable_view()``. Users wishing to “break the chain” will call
``.copy_if_needed()``.

Extended proposal
-----------------

In principle, there’s nothing special about indexing when it comes to defensive
copying. Any method that returns a new ``NDFrame`` without altering existing data
(rename, set_index, possibly assign, dropping columns, etc.) is a candidate for
returning a view. That said, we think it’d be unfortunate if something like the
following was the behavior

.. code-block:: python

   >>> df2 = df.rename(lambda x: x) # suppose df2 is a view on df
   >>> df2.iloc[0, 0] = 10
   MutableViewError("This DataFrame is a view on another DataFrame. Set .as_mutable_view() or copy with ".copy_if_needed()"")

Now we have to ask: does a reasonable consumer of the pandas API expect ``df2``
to be a view? Such that mutating ``df2`` would mutate ``df``? I’d argue no,
people wouldn’t expect that. If that’s the case, then I think requiring people
to include a ``.as_mutable_view()`` or ``.copy_if_needed()`` would be unfortunate
line noise. So in this extended proposal we would probably prefer Copy-on-Write
over Error-on-Write. That said, we don’t wish to discuss the extended proposal
much here. We wish to focus primarily on indexing, and we can make a choice that
is best for indexing. We only mention it here to inform our choice of
Copy-on-Write vs. Error-on-Write.

Propagating mutation forwards
-----------------------------

Thus far we’ve considered the (more common) case of taking a subset, mutating
the subset, and how that should affect the parent. What about the other
direction, where the parent is mutated?

.. code-block:: python

   >>> df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
   >>> df2 = df[["A"]]
   >>> df.iloc[0, 0] = 10
   >>> df2.iloc[0, 0] # what is this value?

We might value symmetry with the “backwards” case, which would argue that the
setitem above should raise (under Error on Write) or copy (under Copy on Write).
Users wishing that setitem operations on the parent should propagate to the
child would need to call .as_mutable_view().

Deprecation or breaking change?
-------------------------------

Because of the subtleties around views vs. copies and mutation, we propose doing
this as an API breaking change accompanying a major version bump. We think that
simply always returning a view is too large a behavior change (even if the
current semantics aren’t well tested / documented, people have written code
that’s tailored to the current implementation). We also think a deprecation
warning is too noisy. Indexing is too common an operation to include a warning
(even if we limit it to just those operations that previously returned copies).

Interaction with BlockManager, ArrayManager, and Consolidation
--------------------------------------------------------------

This proposal is consistent with either the BlockManager or a proposed
ArrayManager. However, there is a subtle interaction with the BlockManager’s
*inplace* consolidation. Today, some operations (e.g. reductions) perform an
inplace consolidation

.. code-block:: python

   >>> df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
   >>> df2 = df1[["A"]].as_mutable_view() # df2 is a view
   >>> df2.mean() # mean consolidates inplace, causing a copy, breaking the view.
   >>> df2.iloc[0, 0] = 1

It would be unfortunate if the presence or absence of a .mean() call changed the
behavior of the later setitem. We likely have the tools to detect these cases
and warn or raise if they occur. But this proposal would likely work better with
a modified BlockManager that doesn’t do inplace consolidation. This will cause
apparent regressions in the performance for workloads that do indexing followed
by many operations that benefit from consolidation. We might consider exposing
consolidation in the public API, though the details of that are left for a
separate discussion.

This proposal is consistent with the proposed ArrayManager.

Background: Current behaviour of views vs copy
----------------------------------------------

To the best of our knowledge, indexing operations currently return views in the
following cases:

Selecting a single column (as a Series) out of a DataFrame is always a view
(``df['a']``) Slicing columns from a DataFrame creating a subset DataFrame
(``df[['a':'b']]`` or ``df.loc[:, 'a': 'b']``) is a view if the the original
DataFrame consists of a single block (single dtype, consolidated) and if you are
slicing (so not a list selection). In all other cases, getting a subset is
always a copy. Slicing rows can return a view, when the row indexer is a slice
object.

Remaining operations (subsetting rows with a list indexer or boolean mask) in
practice return a copy, and we will raise a ``SettingWithCopy`` warning when the
user tries to modify the subset.

Background: Previous attempts
-----------------------------

We’ve discussed this general issue before.
https://github.com/pandas-dev/pandas/issues/10954 and a few pull requests
(https://github.com/pandas-dev/pandas/pull/12036,
https://github.com/pandas-dev/pandas/pull/11207,
https://github.com/pandas-dev/pandas/pull/11500).
