from __future__ import annotations

import math
import random as rnd

import pytest

import dask.bag as db
from dask.bag import random


def test_choices_size_exactly_k():
    seq = range(20)
    sut = db.from_sequence(seq, npartitions=3)
    li = list(random.choices(sut, k=2).compute())
    assert len(li) == 2
    assert all(i in seq for i in li)


def test_choices_k_bigger_than_bag_size():
    seq = range(3)
    sut = db.from_sequence(seq, npartitions=3)
    li = list(random.choices(sut, k=4).compute())
    assert len(li) == 4
    assert all(i in seq for i in li)


def test_choices_empty_partition():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    sut = sut.repartition(3)
    li = list(random.choices(sut, k=2).compute())
    assert sut.map_partitions(len).compute() == (9, 0, 1)
    assert len(li) == 2
    assert all(i in seq for i in li)


def test_choices_k_bigger_than_smallest_partition_size():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    li = list(random.choices(sut, k=2).compute())
    assert sut.map_partitions(len).compute() == (9, 1)
    assert len(li) == 2
    assert all(i in seq for i in li)


def test_choices_k_equal_bag_size_with_unbalanced_partitions():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    li = list(random.choices(sut, k=10).compute())
    assert sut.map_partitions(len).compute() == (9, 1)
    assert len(li) == 10
    assert all(i in seq for i in li)


def test_choices_with_more_bag_partitons():
    # test with npartitions > split_every
    seq = range(100)
    sut = db.from_sequence(seq, npartitions=10)
    li = list(random.choices(sut, k=10, split_every=8).compute())
    assert sut.map_partitions(len).compute() == (10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    assert len(li) == 10
    assert all(i in seq for i in li)


def test_sample_with_more_bag_partitons():
    # test with npartitions > split_every
    seq = range(100)
    sut = db.from_sequence(seq, npartitions=10)
    li = list(random.sample(sut, k=10, split_every=8).compute())
    assert sut.map_partitions(len).compute() == (10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    assert len(li) == 10
    assert all(i in seq for i in li)
    assert len(set(li)) == len(li)


def test_sample_size_exactly_k():
    seq = range(20)
    sut = db.from_sequence(seq, npartitions=3)
    li = list(random.sample(sut, k=2).compute())
    assert sut.map_partitions(len).compute() == (7, 7, 6)
    assert len(li) == 2
    assert all(i in seq for i in li)
    assert len(set(li)) == len(li)


def test_sample_k_bigger_than_bag_size():
    seq = range(3)
    sut = db.from_sequence(seq, npartitions=3)
    # should raise: Sample larger than population or is negative
    with pytest.raises(ValueError, match="Sample larger than population"):
        random.sample(sut, k=4).compute()


def test_sample_empty_partition():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    sut = sut.repartition(3)
    li = list(random.sample(sut, k=2).compute())
    assert sut.map_partitions(len).compute() == (9, 0, 1)
    assert len(li) == 2
    assert all(i in seq for i in li)
    assert len(set(li)) == len(li)


def test_sample_size_k_bigger_than_smallest_partition_size():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    li = list(random.sample(sut, k=2).compute())
    assert sut.map_partitions(len).compute() == (9, 1)
    assert len(li) == 2
    assert all(i in seq for i in li)
    assert len(set(li)) == len(li)


def test_sample_k_equal_bag_size_with_unbalanced_partitions():
    seq = range(10)
    sut = db.from_sequence(seq, partition_size=9)
    li = list(random.sample(sut, k=10).compute())
    assert sut.map_partitions(len).compute() == (9, 1)
    assert len(li) == 10
    assert all(i in seq for i in li)
    assert len(set(li)) == len(li)


def test_sample_k_larger_than_partitions():
    bag = db.from_sequence(range(10), partition_size=3)
    bag2 = random.sample(bag, k=8, split_every=2)
    seq = bag2.compute()
    assert len(seq) == 8


def test_weighted_sampling_without_replacement():
    population = range(4)
    p = [0.01, 0.33, 0.33, 0.33]
    k = 3
    sampled = random._weighted_sampling_without_replacement(
        population=population, weights=p, k=k
    )
    assert len(set(sampled)) == k


def test_sample_return_bag():
    seq = range(20)
    sut = db.from_sequence(seq, npartitions=3)
    assert isinstance(random.sample(sut, k=2), db.Bag)


def test_partitions_are_coerced_to_lists():
    # https://github.com/dask/dask/issues/6906
    A = db.from_sequence([[1, 2], [3, 4, 5], [6], [7]])
    B = db.from_sequence(["a", "b", "c", "d"])

    a = random.choices(A.flatten(), k=B.count().compute()).repartition(4)

    C = db.zip(B, a).compute()
    assert len(C) == 4


def test_reservoir_sample_map_partitions_correctness():
    N, k = 20, 10
    seq = list(range(N))
    distribution = [0 for _ in range(N)]
    expected_distribution = [0 for _ in range(N)]
    reps = 2000
    for _ in range(reps):
        picks, _ = random._sample_map_partitions(seq, k)
        for pick in picks:
            distribution[pick] += 1
        for pick in rnd.sample(seq, k=k):
            expected_distribution[pick] += 1

    # convert to probabilities
    distribution = [c / (reps * k) for c in distribution]
    expected_distribution = [c / (reps * k) for c in expected_distribution]

    # use bhattacharyya distance to asses the similarity of distributions
    assert math.isclose(
        0.0, bhattacharyya(distribution, expected_distribution), abs_tol=1e-2
    )


def test_reservoir_sample_with_replacement_map_partitions_correctness():
    N, k = 20, 10
    seq = list(range(N))
    distribution = [0 for _ in range(N)]
    expected_distribution = [0 for _ in range(N)]
    reps = 2000
    for _ in range(reps):
        picks, _ = random._sample_with_replacement_map_partitions(seq, k)
        for pick in picks:
            distribution[pick] += 1
        for pick in rnd.choices(seq, k=k):
            expected_distribution[pick] += 1

    # convert to probabilities
    distribution = [c / (reps * k) for c in distribution]
    expected_distribution = [c / (reps * k) for c in expected_distribution]

    # use bhattacharyya distance to asses the similarity of distributions
    assert math.isclose(
        0.0, bhattacharyya(distribution, expected_distribution), abs_tol=1e-2
    )


def bhattacharyya(h1, h2):
    return 1 - sum(math.sqrt(hi * hj) for hi, hj in zip(h1, h2))
