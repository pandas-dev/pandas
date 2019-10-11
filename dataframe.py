import pandas as pd
import numpy as np
import sys

data = pd.DataFrame(np.array([[1, 2, 3], [4, np.NaN, 6], [7, 8, 9]]),
                    columns=['asdfadsffdsafdsafsaa', 'b', 'c'])
print(data.info(verbose=True))
print(data)
print(data.info())