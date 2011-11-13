from pandas import *
import string


n = 100000
asize = 5
bsize = 5

letters = np.asarray(list(string.letters), dtype=object)

data = DataFrame(dict(foo=letters[:asize][np.random.randint(0, asize, n)],
                      bar=letters[:bsize][np.random.randint(0, bsize, n)],
                      baz=np.random.randn(n),
                      qux=np.random.randn(n)))

table = pivot_table(data, xby=['foo', 'bar'])
