import numpy as np
import pandas
from pandas import Series,DataFrame

print pandas.__version__

s = Series(np.arange(1028.))

df = DataFrame({ i:s for i in range(1028) })

import pdb; pdb.set_trace()
df.apply(lambda x: np.corrcoef(x,s)[0,1])
