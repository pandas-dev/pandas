import pandas as pd
from datetime import datetime
df = pd.Series(
    [datetime(2000, 1, 1)], dtype="timestamp[ns][pyarrow]", name="col0"
).to_frame()

dfi = df.__dataframe__()
result = pd.api.interchange.from_dataframe(dfi)

print(result)
# 1970-01-02 14:55:26.266485584