"""
Interaction with scipy.sparse matrices.

Currently only includes SparseSeries.to_coo helpers.
"""
from pandas.core.frame import DataFrame
from pandas.core.index import MultiIndex, Index
from pandas.core.series import Series
import itertools
import numpy
from pandas.compat import OrderedDict
from pandas.tools.util import cartesian_product


def _get_label_to_i_dict(labels, sort_labels=False):
    """ Return OrderedDict of unique labels to number. Optionally sort by label. """
    labels = Index(map(tuple, labels)).unique().tolist()  # squish
    if sort_labels:
        labels = sorted(list(labels))
    d = OrderedDict((k, i) for i, k in enumerate(labels))
    return(d)


def _get_index_subset_to_coord_dict(index, subset, sort_labels=False):
    ilabels = list(zip(*[index.get_level_values(i) for i in subset]))
    labels_to_i = _get_label_to_i_dict(ilabels, sort_labels=sort_labels)
    return(labels_to_i)


def _check_is_partition(parts, whole):
    whole = set(whole)
    parts = [set(x) for x in parts]
    if set.intersection(*parts) != set():
        raise ValueError(
            'Is not a partition because intersection is not null.')
    if set.union(*parts) != whole:
        raise ValueError('Is not a partition becuase union is not the whole.')


def _to_ijv(ss, ilevels=(0,), jlevels=(1,), sort_labels=False):
    """ For arbitrary (MultiIndexed) SparseSeries return
    (v, i, j, ilabels, jlabels) where (v, (i, j)) is suitable for
    passing to scipy.sparse.coo constructor. """
    # index and column levels must be a partition of the index
    _check_is_partition([ilevels, jlevels], range(ss.index.nlevels))

    # from the SparseSeries: get the labels and data for non-null entries
    values = ss._data.values._valid_sp_values
    blocs = ss._data.values.sp_index.blocs
    blength = ss._data.values.sp_index.blengths
    nonnull_labels = list(
        itertools.chain(*[ss.index.values[i:(i + j)] for i, j in zip(blocs, blength)]))

    def get_indexers(levels):
        """ Return sparse coords and dense labels for subset levels """
        values_ilabels = [tuple(x[i] for i in levels) for x in nonnull_labels]
        labels_to_i = _get_index_subset_to_coord_dict(
            ss.index, levels, sort_labels=sort_labels)
        i_coord = [labels_to_i[i] for i in values_ilabels]
        return(i_coord, list(labels_to_i.keys()))

    i_coord, i_labels = get_indexers(ilevels)
    j_coord, j_labels = get_indexers(jlevels)

    return(values, i_coord, j_coord, i_labels, j_labels)


def _sparse_series_to_coo(ss, ilevels=(0,), jlevels=(1,), sort_labels=False):
    """ Convert a SparseSeries to a scipy.sparse.coo_matrix using index levels ilevels, jlevels as the row and column
    labels respectively. Returns the sparse_matrix, row and column labels. """
    if ss.index.nlevels < 2:
        raise ValueError('to_coo requires MultiIndex with nlevels > 2')
    if not ss.index.is_unique:
        raise ValueError(
            'Duplicate index entries are not allowed in to_coo transformation.')

    # to keep things simple, only rely on integer indexing (not labels)
    ilevels = [ss.index._get_level_number(x) for x in ilevels]
    jlevels = [ss.index._get_level_number(x) for x in jlevels]
    ss = ss.copy()
    ss.index.names = [None] * ss.index.nlevels  # kill any existing labels

    v, i, j, il, jl = _to_ijv(
        ss, ilevels=ilevels, jlevels=jlevels, sort_labels=sort_labels)
    import scipy.sparse
    sparse_matrix = scipy.sparse.coo_matrix(
        (v, (i, j)), shape=(len(il), len(jl)))
    return(sparse_matrix, il, jl)


def _coo_to_sparse_series(A, dense_index=False):
    """ Convert a scipy.sparse.coo_matrix to a SparseSeries.
    Use the defaults given in the SparseSeries constructor. """
    s = Series(A.data, MultiIndex.from_arrays((A.row, A.col)))
    s = s.sort_index()
    s = s.to_sparse()  # TODO: specify kind?
    if dense_index:
        # is there a better constructor method to use here?
        i = range(A.shape[0])
        j = range(A.shape[1])
        ind = MultiIndex.from_product([i, j])
        s = s.reindex_axis(ind)
    return(s)
