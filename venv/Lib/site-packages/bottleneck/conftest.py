
import pytest

try:
    import pytest_run_parallel  # noqa:F401

    PARALLEL_RUN_AVAILABLE = True
except Exception:
    PARALLEL_RUN_AVAILABLE = False


def pytest_configure(config):
    if not PARALLEL_RUN_AVAILABLE:
        config.addinivalue_line(
            'markers',
            'parallel_threads(n): run the given test function in parallel '
            'using `n` threads.',
        )
        config.addinivalue_line(
            "markers",
            "thread_unsafe: mark the test function as single-threaded",
        )
        config.addinivalue_line(
            "markers",
            "iterations(n): run the given test function `n` times in each thread",
        )


if not PARALLEL_RUN_AVAILABLE:
    @pytest.fixture
    def num_parallel_threads():
        return 1
