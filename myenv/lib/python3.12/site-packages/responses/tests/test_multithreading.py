"""
Separate file for multithreading since it takes time to run
"""
import threading

import pytest
import requests

import responses


@pytest.mark.parametrize("execution_number", range(10))
def test_multithreading_lock(execution_number):  # type: ignore[misc]
    """Reruns test multiple times since error is random and
    depends on CPU and can lead to false positive result.

    """
    n_threads = 10
    n_requests = 30
    with responses.RequestsMock() as m:
        for j in range(n_threads):
            for i in range(n_requests):
                m.add(url=f"http://example.com/example{i}", method="GET")

        def fun():
            for req in range(n_requests):
                requests.get(f"http://example.com/example{req}")

        threads = [
            threading.Thread(name=f"example{i}", target=fun) for i in range(n_threads)
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
