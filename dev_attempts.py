import pandas as pd
import numpy as np
import timeit
np.random.seed(43)
# also tested for n = 1000, 10_000, 100_000
n=1_000_000
cols = list('abcdef')
df = pd.DataFrame(np.random.randint(0, 10, size=(n,len(cols))), columns=cols)
df['col'] = np.random.choice(cols, n)
idx = df['col'].index.to_numpy()
cols = df['col'].to_numpy()

def og_lookup(idx, cols):
	return df.lookup(idx, cols,'og')

# def melt_lookup():
# 	melt = df.melt('col')
# 	melt = melt.loc[lambda x: x['col']==x['variable'], 'value']
# 	melt = melt.reset_index(drop=True)
# 	return melt

# def quan_lookup(idx,cols):
# 	return df.reindex(cols,axis=1).to_numpy()[np.arange(df.shape[0]), idx]

# def quan_lookup2(idx,cols):
# 	return df.reindex(cols,axis=1).to_numpy()[np.arange(df.shape[0]), idx]

# def marco_lookup():
# 	return df.melt('col', ignore_index=False).query('col==variable')['value'].reindex(df.index).to_numpy()


timeit.timeit(lambda: og_lookup(idx,cols),number=10)
# timeit.timeit(lambda: melt_lookup(idx,cols),number=10)
# timeit.timeit(lambda: quan_lookup(idx,cols),number=10)
# timeit.timeit(lambda: quan_lookup2(idx,cols),number=10)
# timeit.timeit(lambda: marco_lookup(idx,cols),number=10)

# idx, cols = pd.factorize(df['col'])
# df.reindex(cols, axis=1).to_numpy()[np.arange(len(df)), idx]
