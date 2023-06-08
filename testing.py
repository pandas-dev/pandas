import numpy as np

import pandas as pd

comment_lines = ["Hello", "world", "longer sentence here"]

df = pd.DataFrame(np.random.rand(4, 4))

df.to_csv("testing.csv", index=False)
df.to_csv("testing.csv2", index=False, comment="#", comment_lines=comment_lines)
