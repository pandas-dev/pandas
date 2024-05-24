import pandas as pd
import multiprocessing as mp

def parallel_apply(df, func, num_partitions=None):
    """
    Apply a function to a DataFrame in parallel using multiprocessing.
    
    Parameters:
    - df: Pandas DataFrame
    - func: Function to apply to each partition
    - num_partitions: Number of partitions to split the DataFrame into (default is the number of CPUs)
    """
    if num_partitions is None:
        num_partitions = mp.cpu_count()
    
    df_split = np.array_split(df, num_partitions)
    pool = mp.Pool(num_partitions)
    result = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    
    return result

data = {'A': range(1000000), 'B': range(1000000)}
df = pd.DataFrame(data)
# result = parallel_apply(df, lambda x: x**2)