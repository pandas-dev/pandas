"""
K-means clustering and vector quantization (:mod:`scipy.cluster.vq`)
====================================================================

Provides routines for k-means clustering, generating code books
from k-means models and quantizing vectors by comparing them with
centroids in a code book.

.. autosummary::
   :toctree: generated/

   whiten -- Normalize a group of observations so each feature has unit variance
   vq -- Calculate code book membership of a set of observation vectors
   kmeans -- Perform k-means on a set of observation vectors forming k clusters
   kmeans2 -- A different implementation of k-means with more methods
           -- for initializing centroids

Exceptions:

.. autosummary::
   :toctree: generated/

   ClusterError

Background information
----------------------
The k-means algorithm takes as input the number of clusters to
generate, k, and a set of observation vectors to cluster. It
returns a set of centroids, one for each of the k clusters. An
observation vector is classified with the cluster number or
centroid index of the centroid closest to it.

A vector v belongs to cluster i if it is closer to centroid i than
any other centroid. If v belongs to i, we say centroid i is the
dominating centroid of v. The k-means algorithm tries to
minimize distortion, which is defined as the sum of the squared distances
between each observation vector and its dominating centroid.
The minimization is achieved by iteratively reclassifying
the observations into clusters and recalculating the centroids until
a configuration is reached in which the centroids are stable. One can
also define a maximum number of iterations.

Since vector quantization is a natural application for k-means,
information theory terminology is often used. The centroid index
or cluster index is also referred to as a "code" and the table
mapping codes to centroids and, vice versa, is often referred to as a
"code book". The result of k-means, a set of centroids, can be
used to quantize vectors. Quantization aims to find an encoding of
vectors that reduces the expected distortion.

All routines expect obs to be an M by N array, where the rows are
the observation vectors. The codebook is a k by N array, where the
ith row is the centroid of code word i. The observation vectors
and centroids have the same feature dimension.

As an example, suppose we wish to compress a 24-bit color image
(each pixel is represented by one byte for red, one for blue, and
one for green) before sending it over the web. By using a smaller
8-bit encoding, we can reduce the amount of data by two
thirds. Ideally, the colors for each of the 256 possible 8-bit
encoding values should be chosen to minimize distortion of the
color. Running k-means with k=256 generates a code book of 256
codes, which fills up all possible 8-bit sequences. Instead of
sending a 3-byte value for each pixel, the 8-bit centroid index
(or code word) of the dominating centroid is transmitted. The code
book is also sent over the wire so each 8-bit code can be
translated back to a 24-bit pixel value representation. If the
image of interest was of an ocean, we would expect many 24-bit
blues to be represented by 8-bit codes. If it was an image of a
human face, more flesh-tone colors would be represented in the
code book.

"""
from ._vq_impl import ClusterError, kmeans, kmeans2, vq, whiten
from ._vq_impl import py_vq

__all__ = ["ClusterError", "kmeans", "kmeans2", "vq", "whiten"]

# deprecated attributes
__all__ += ["py_vq"]


def __dir__(): 
    return __all__ 
