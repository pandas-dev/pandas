df = pd.DataFrame({'lab': ['A', 'B', 'C'], 'val': [10, 30, 20]})
ax = df.plot.barh(x='lab', y='val')
