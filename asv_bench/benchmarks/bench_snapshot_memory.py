import pandas as pd
import numpy as np
import tracemalloc

def make_df(nrows=1_000_000):
    return pd.DataFrame({
        "a": np.random.randint(0, 100, size=nrows),
        "b": np.random.random(size=nrows),
        "c": np.random.choice(list("abcdefghijklmnopqrstuvwxyz"), size=nrows)
    })

df = make_df(200_000)
tracemalloc.start()
snap = df.snapshot("bench")
snap_shot = tracemalloc.take_snapshot()
top_stats = snap_shot.statistics('lineno')
print("Top memory stats:", top_stats[:3])
