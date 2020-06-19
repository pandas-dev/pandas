axes = df.plot.bar(
    rot=0, subplots=True, color={"speed": "red", "lifespan": "green"}
)
axes[1].legend(loc=2)  # doctest: +SKIP
