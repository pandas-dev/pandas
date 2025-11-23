import io

import pandas as pd

print("Testing issue #63089...")
try:
    result = pd.read_csv(
        io.StringIO("""h
4e492493924924""")
    )
    print("Success! Result:")
    print(result)
except Exception as e:
    print(f"Exception occurred: {type(e).__name__}: {e}")
