df = pd.DataFrame(np.random.randn(10, 3),
                  columns=['Col1', 'Col2', 'Col3'])
df['X'] = pd.Series(['A', 'A', 'A', 'A', 'A',
                     'B', 'B', 'B', 'B', 'B'])
df['Y'] = pd.Series(['A', 'B', 'A', 'B', 'A',
                     'B', 'A', 'B', 'A', 'B'])
boxplot = df.boxplot(column=['Col1', 'Col2'], by=['X', 'Y'])
