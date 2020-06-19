n = 500
df = pd.DataFrame({
    'coord_x': np.random.uniform(-3, 3, size=n),
    'coord_y': np.random.uniform(30, 50, size=n),
    'observations': np.random.randint(1,5, size=n)
    })
ax = df.plot.hexbin(x='coord_x',
                    y='coord_y',
                    C='observations',
                    reduce_C_function=np.sum,
                    gridsize=10,
                    cmap="viridis")
