# 0 - 4 parameters (`'gausshyper'`)
distcont: list[
    tuple[str, tuple[()] | tuple[float] | tuple[float, float] | tuple[float, float, float] | tuple[float, float, float, float]]
]

# 0 - 4 parameters (`'gausshyper'`)
invdistcont: list[
    tuple[str, tuple[()] | tuple[float] | tuple[float, float] | tuple[float, float, float] | tuple[float, float, float, float]]
]

# 1 - 4 parameters (`'nchypergeom_fisher'` and `'nchypergeom_wallenius'`)
distdiscrete: list[tuple[str, tuple[float] | tuple[float, float] | tuple[int, float, float] | tuple[int, int, int, float]]]
# 1 - 4 parameters (`'nchypergeom_fisher'` and `'nchypergeom_wallenius'`)
invdistdiscrete: list[tuple[str, tuple[float] | tuple[float, float] | tuple[int, float, float] | tuple[int, int, int, float]]]
