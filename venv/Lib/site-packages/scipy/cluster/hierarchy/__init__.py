"""
Hierarchical clustering (:mod:`scipy.cluster.hierarchy`)
========================================================

.. currentmodule:: scipy.cluster.hierarchy

These functions cut hierarchical clusterings into flat clusterings
or find the roots of the forest formed by a cut by providing the flat
cluster ids of each observation.

.. autosummary::
   :toctree: generated/

   fcluster
   fclusterdata
   leaders

These are routines for agglomerative clustering.

.. autosummary::
   :toctree: generated/

   linkage
   single
   complete
   average
   weighted
   centroid
   median
   ward

These routines compute statistics on hierarchies.

.. autosummary::
   :toctree: generated/

   cophenet
   from_mlab_linkage
   inconsistent
   maxinconsts
   maxdists
   maxRstat
   to_mlab_linkage

Routines for visualizing flat clusters.

.. autosummary::
   :toctree: generated/

   dendrogram

These are data structures and routines for representing hierarchies as
tree objects.

.. autosummary::
   :toctree: generated/

   ClusterNode
   leaves_list
   to_tree
   cut_tree
   optimal_leaf_ordering

These are predicates for checking the validity of linkage and
inconsistency matrices as well as for checking isomorphism of two
flat cluster assignments.

.. autosummary::
   :toctree: generated/

   is_valid_im
   is_valid_linkage
   is_isomorphic
   is_monotonic
   correspond
   num_obs_linkage

Utility routines for plotting:

.. autosummary::
   :toctree: generated/

   set_link_color_palette

Utility classes:

.. autosummary::
   :toctree: generated/

   DisjointSet -- data structure for incremental connectivity queries
   
Warnings:
    
.. autosummary::
   :toctree: generated/

   ClusterWarning

"""
from ._hierarchy_impl import (
    ClusterNode, ClusterWarning, DisjointSet, average, centroid, complete, cophenet,
    correspond, cut_tree, dendrogram, fcluster, fclusterdata, from_mlab_linkage,
    inconsistent, is_isomorphic, is_monotonic, is_valid_im, is_valid_linkage, leaders,
    leaves_list, linkage, maxRstat, maxdists, maxinconsts, median, num_obs_linkage,
    optimal_leaf_ordering, set_link_color_palette, single, to_mlab_linkage, to_tree,
    ward, weighted
)

__all__ = ['ClusterNode', 'ClusterWarning', 'DisjointSet', 'average', 'centroid',
           'complete', 'cophenet', 'correspond', 'cut_tree', 'dendrogram', 'fcluster',
           'fclusterdata', 'from_mlab_linkage', 'inconsistent',
           'is_isomorphic', 'is_monotonic', 'is_valid_im', 'is_valid_linkage',
           'leaders', 'leaves_list', 'linkage', 'maxRstat', 'maxdists',
           'maxinconsts', 'median', 'num_obs_linkage', 'optimal_leaf_ordering',
           'set_link_color_palette', 'single', 'to_mlab_linkage', 'to_tree',
           'ward', 'weighted']
