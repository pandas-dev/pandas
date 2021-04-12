import pytest
import pandas as pd
from datetime import datetime as dt

"""targeting #GH40781"""


def test_if_plotable_xlim_ylim_both_ints() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[1, 2], ylim=[3, 4])
        return True
    except ValueError:
        pytest.raises(ValueError, match="FAILED")


def test_if_plotable_xlim_first_int() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[1, "2"], ylim=[3, 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_xlim_first_str() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=["1", 2], ylim=[3, 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_ylim_first_int() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[1, 2], ylim=[3, "4"])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_ylim_first_str() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[1, 2], ylim=["3", 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_ylim_both_str() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[1, 2], ylim=["3", "4"])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_xlim_both_str() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=["1", "2"], ylim=[3, 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_xlim_datetime() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[dt.now(), dt.now()], ylim=[3, 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True


def test_if_plotable_xlim_datetime_str() -> bool:
    # checks if ValueError is raised with invalid dtype
    df = pd.DataFrame({
        "A": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "B" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    })
    try:
        df.plot(kind="line", x="A", xlim=[dt.now(), "1"], ylim=[3, 4])
        pytest.raises(ValueError, match="FAILED")
    except ValueError:
        return True
