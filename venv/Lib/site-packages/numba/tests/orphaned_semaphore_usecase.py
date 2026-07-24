# See issue: https://github.com/numba/numba/issues/4348
# this file must be run in isolation to replicate, test:
# numba.tests.test_parallel_backend.TestInitSafetyIssues.test_orphaned_semaphore
# does this to check semaphores are not leaking.

import multiprocessing as mp

import numba # noqa, deliberately unused, here to test import is safe


def w():
    pass


def main():
    ps = [mp.Process(target=w) for _ in range(4)]
    [p.start() for p in ps]
    [p.join() for p in ps]


if __name__ == '__main__':
    p = mp.get_context('spawn').Process(target=main)
    p.start()
    p.join()
