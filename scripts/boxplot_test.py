import matplotlib.pyplot as plt

import random
import pandas.util.testing as tm
tm.N = 1000
df = tm.makeTimeDataFrame()
import string
foo = list(string.letters[:5]) * 200
df['indic'] = list(string.letters[:5]) * 200
random.shuffle(foo)
df['indic2'] = foo
df.boxplot(by=['indic', 'indic2'], fontsize=8, rot=90)

plt.show()
