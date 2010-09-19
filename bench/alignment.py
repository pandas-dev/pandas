# Setup
import numpy as np
import pandas
import la
N = 1000
K = 50
arr1 = np.random.randn(N, K)
arr2 = np.random.randn(N, K)
idx1 = range(N)
idx2 = range(K)

# pandas
dma1 = pandas.DataMatrix(arr1, idx1, idx2)
dma2 = pandas.DataMatrix(arr2, idx1[::-1], idx2[::-1])

# larry
lar1 = la.larry(arr1, [idx1, idx2])
lar2 = la.larry(arr2, [idx1[::-1], idx2[::-1]])

for i in range(100):
    result = lar1 + lar2
