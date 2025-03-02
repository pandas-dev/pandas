import numpy as np

import pandas as pd

# creating a single series dataframe
frame = pd.DataFrame(np.array([1, 2, 3, 4, 5, 100]))

# getting the describe with single percentile value
print(frame.describe(percentiles=[0.25]))
