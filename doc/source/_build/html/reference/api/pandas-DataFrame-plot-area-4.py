df = pd.DataFrame({
    'sales': [3, 2, 3],
    'visits': [20, 42, 28],
    'day': [1, 2, 3],
})
ax = df.plot.area(x='day')
