import pandas as pd

a = 1450804465901089690
b = 1450804465901089614

aS = pd.Series([a])
bS = pd.Series([b])

atol = 100
rtol = 0

diff = abs(a-b)
if diff > atol + rtol*abs(b): 
    # Formula came from https://pandas.pydata.org/docs/reference/api/pandas.testing.assert_frame_equal.html > check_exact argument description
    print("Substantial difference detected:")
    print(diff)
else:
    print("No substantial difference was detected...and yet...")
    

pd.testing.assert_series_equal(aS, bS, check_dtype=False, check_exact=False, atol=atol, rtol=rtol)