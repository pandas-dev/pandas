import pandas as pd
import numpy as np
import timeit
np.random.seed(43)
for n in [100,100_000]:
	cols = list('abcdef')
	df = pd.DataFrame(np.random.randint(0, 10, size=(n,len(cols))), columns=cols)
	df['col'] = np.random.choice(cols, n)
	idx = df['col'].index.to_numpy()
	cols = df['col'].to_numpy()
	timeit.timeit(lambda: df.lookup(idx, cols,'og'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'a'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'b'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'c'),number=10)
	df['a'] = df['a'].astype(str)
	df['a'] = 'a'
	print('mixed')
	timeit.timeit(lambda: df.lookup(idx, cols,'og'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'a'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'b'),number=10)
	timeit.timeit(lambda: df.lookup(idx, cols,'c'),number=10)
	df.lookup(idx, cols,'b')
	print('\n')