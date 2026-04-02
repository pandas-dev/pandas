import pandas as pd
import numpy as np
import timeit
from pandas.core.internals.managers import BlockManager

def benchmark_is_consolidated():
    # Create a DataFrame with many columns of different dtypes to trigger 
    # many blocks in the BlockManager.
    # 50 float columns, 50 int columns -> 100 blocks total (pre-consolidation)
    data = {}
    for i in range(500):
        data[f'float_{i}'] = np.random.randn(100)
        data[f'int_{i}'] = np.random.randint(0, 100, 100)
    
    df = pd.DataFrame(data)
    
    # Access the internal BlockManager
    mgr = df._mgr
    
    from pandas.core.internals.blocks import new_block
    from pandas.core.internals.managers import BlockManager
    from pandas.core.internals.construction import BlockPlacement
    
    blocks = []
    for i, (col, vals) in enumerate(data.items()):
        bp = BlockPlacement(slice(i, i+1))
        blk = new_block(vals.reshape(1, -1), placement=bp, ndim=2)
        blocks.append(blk)
    
    mgr = BlockManager(tuple(blocks), [pd.Index(data.keys()), pd.Index(range(100))])
    
    def run_check():
        mgr._known_consolidated = False
        return mgr.is_consolidated()

    # Benchmark the execution time
    timer = timeit.Timer(run_check)
    iterations = 10000
    total_time = timer.timeit(number=iterations)
    
    print(f"Total time for {iterations} calls to is_consolidated(): {total_time:.4f} seconds")
    print(f"Average time per call: {total_time/iterations*1e6:.4f} microseconds")

if __name__ == "__main__":
    benchmark_is_consolidated()
