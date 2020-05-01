import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.DataFrame(dict(x=[np.nan, np.nan]))
    import pdb; pdb.set_trace()
    res = df.max()
    print(res)
