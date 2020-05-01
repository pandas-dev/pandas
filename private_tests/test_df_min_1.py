import pandas as pd

if __name__ == '__main__':
    df = pd.DataFrame(dict(x=pd.to_datetime([])))
    series = pd.Series(pd.to_datetime([]))
    res = df.min(1)
    print(res)

