df = pd.read_csv(
    'https://raw.github.com/pandas-dev/'
    'pandas/master/pandas/tests/data/iris.csv'
)
pd.plotting.andrews_curves(df, 'Name')
