# Let's repeat the same example, but specifying colors for
# each column (in this case, for each animal).

axes = df.plot.line(
    subplots=True, color={"pig": "pink", "horse": "#742802"}
)
