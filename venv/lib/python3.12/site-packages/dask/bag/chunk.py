from __future__ import annotations


def barrier(*args):
    return None


def getitem(x, key):
    """Like :func:`operator.getitem`, but allows setting key using partial
    ``partial(chunk.getitem, key=key)
    """
    return x[key]


def foldby_combine2(combine, acc, x):
    return combine(acc, x[1])


def groupby_tasks_group_hash(x, hash, grouper):
    return hash(grouper(x)), x


def var_chunk(seq):
    squares, total, n = 0.0, 0.0, 0
    for x in seq:
        squares += x**2
        total += x
        n += 1
    return squares, total, n


def var_aggregate(x, ddof):
    squares, totals, counts = list(zip(*x))
    x2, x, n = float(sum(squares)), float(sum(totals)), sum(counts)
    result = (x2 / n) - (x / n) ** 2
    return result * n / (n - ddof)
