# Setup
from pandas.compat import range, lrange
import numpy as np
import pandas
import la
N = 1000
K = 50
arr1 = np.random.randn(N, K)
arr2 = np.random.randn(N, K)
idx1 = lrange(N)
idx2 = lrange(K)

# pandas
dma1 = pandas.DataFrame(arr1, idx1, idx2)
dma2 = pandas.DataFrame(arr2, idx1[::-1], idx2[::-1])

# larry
lar1 = la.larry(arr1, [idx1, idx2])
lar2 = la.larry(arr2, [idx1[::-1], idx2[::-1]])

for i in range(100):
    result = lar1 + lar2
