import pandas as pd

def test_empty_resample():
    # Test 1: Empty DataFrame with datetime index
    df1 = pd.DataFrame(index=pd.date_range('2024-01-01', '2024-01-02', freq='h'))
    print("\nTest 1: Empty DataFrame with datetime index")
    print("Original index type:", type(df1.index))
    resampled1 = df1.resample('2h').mean()
    print("Resampled index type:", type(resampled1.index))
    
    # Test 2: Empty DataFrame with multi-index
    df2 = pd.DataFrame(index=pd.MultiIndex.from_product([
        pd.date_range('2024-01-01', '2024-01-02', freq='h'),
        ['A', 'B']
    ]))
    print("\nTest 2: Empty DataFrame with multi-index")
    print("Original index type:", type(df2.index))
    resampled2 = df2.resample('2h', level=0).mean()
    print("Resampled index type:", type(resampled2.index))
    
    # Test 3: Empty DataFrame with columns but no data
    df3 = pd.DataFrame(columns=['A', 'B'], index=pd.date_range('2024-01-01', '2024-01-02', freq='h'))
    print("\nTest 3: Empty DataFrame with columns but no data")
    print("Original index type:", type(df3.index))
    resampled3 = df3.resample('2h').mean()
    print("Resampled index type:", type(resampled3.index))

if __name__ == "__main__":
    test_empty_resample() 