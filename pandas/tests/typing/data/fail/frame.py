# flake8: noqa

import pandas as pd

empty_df = pd.DataFrame()
empty_ser = pd.Series()

empty_df.dot()  # E: All overload variants of "dot" of "DataFrame" require at least one argument
