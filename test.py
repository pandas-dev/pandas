

import pandas as pd
df = pd.DataFrame(
    {'pid' : [1,1,1,2,2,3,3,3],
     'tag' : [23,45,62,24,45,34,25,62],
          })

g = df.groupby('tag')

import pdb; pdb.set_trace()
g.filter(lambda x: len(x) > 1)
