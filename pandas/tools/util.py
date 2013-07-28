from pandas.core.index import Index
import numpy as np

def match(needles, haystack):
    haystack = Index(haystack)
    needles = Index(needles)
    return haystack.get_indexer(needles)

def cartesian_product(X):
    '''
    Numpy version of itertools.product or pandas.compat.product.
    Sometimes faster (for large inputs)...

    Examples
    --------
    >>> cartesian_product([list('ABC'), [1, 2]])
    [array(['A', 'A', 'B', 'B', 'C', 'C'], dtype='|S1'),
 	array([1, 2, 1, 2, 1, 2])]

    '''

    lenX = np.fromiter((len(x) for x in X), dtype=int)
    cumprodX = np.cumproduct(lenX)

    a = np.roll(cumprodX, 1)
    a[0] = 1

    b = cumprodX[-1] / cumprodX

    return [np.tile(np.repeat(x, b[i]), 
                    np.product(a[i]))
               for i, x in enumerate(X)]