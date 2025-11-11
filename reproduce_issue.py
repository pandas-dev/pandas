import numpy as np
import pandas as pd
import joblib

def job():
    df = pd.DataFrame({
        'date' : [np.datetime64("20110101", 'ns')], 
        'id' : [1], 
        'val' : [1]
    }).set_index(["date", "id"])
    return df

# This should trigger the error
try:
    result = joblib.Parallel(n_jobs=2)(
        [joblib.delayed(job)()]
    )
    print("No error occurred - SUCCESS!")
    print(result)
except Exception as e:
    print(f"Error occurred: {e}")
    import traceback
    traceback.print_exc()