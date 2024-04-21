import pandas as pd
import datetime

# Attempt to replicate the issue
result = pd.DataFrame([datetime.datetime.now(datetime.timezone.utc)])[0].value_counts(bins=2)
print(result)
print(pd.__version__)
