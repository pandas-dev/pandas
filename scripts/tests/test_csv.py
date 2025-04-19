import tempfile

import numpy as np

import pandas as pd


def test_preserve_numpy_arrays_in_csv():
    print("\nRunning: test_preserve_numpy_arrays_in_csv")
    df = pd.DataFrame({
        "id": [1, 2],
        "embedding": [
            np.array([0.1, 0.2, 0.3]),
            np.array([0.4, 0.5, 0.6]),
        ],
    })

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmp:
        path = tmp.name
        df.to_csv(path, index=False, preserve_complex=True)
        df_loaded = pd.read_csv(path, preserve_complex=True)

        assert isinstance(
            df_loaded["embedding"][0], np.ndarray
        ), "Test Failed: The CSV did not preserve embeddings as NumPy arrays!"

    print("PASS: test_preserve_numpy_arrays_in_csv")


def test_preserve_numpy_arrays_in_csv_empty_dataframe():
    print("\nRunning: test_preserve_numpy_arrays_in_csv_empty_dataframe")
    df = pd.DataFrame({"embedding": []})
    expected = "embedding\n"

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmp:
        path = tmp.name
        df.to_csv(path, index=False, preserve_complex=True)
        with open(path, encoding="utf-8") as f:
            result = f.read()

    msg = (
        f"CSV output mismatch for empty DataFrame.\n"
        f"Got:\n{result}\nExpected:\n{expected}"
    )
    assert result == expected, msg
    print("PASS: test_preserve_numpy_arrays_in_csv_empty_dataframe")


def test_preserve_numpy_arrays_in_csv_mixed_dtypes():
    print("\nRunning: test_preserve_numpy_arrays_in_csv_mixed_dtypes")
    df = pd.DataFrame({
        "id": [101, 102],
        "name": ["alice", "bob"],
        "scores": [
            np.array([95.5, 88.0]),
            np.array([76.0, 90.5]),
        ],
        "age": [25, 30],
    })

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmp:
        path = tmp.name
        df.to_csv(path, index=False, preserve_complex=True)
        df_loaded = pd.read_csv(path, preserve_complex=True)

        err_scores = "Failed: 'scores' column not deserialized as np.ndarray."
        assert isinstance(df_loaded["scores"][0], np.ndarray), err_scores
        assert df_loaded["id"].dtype == np.int64, (
            "Failed: 'id' should still be int."
        )
        assert df_loaded["name"].dtype == object, (
            "Failed: 'name' should still be object/string."
        )
        assert df_loaded["age"].dtype == np.int64, (
            "Failed: 'age' should still be int."
        )

    print("PASS: test_preserve_numpy_arrays_in_csv_mixed_dtypes")


if __name__ == "__main__":
    test_preserve_numpy_arrays_in_csv()
    test_preserve_numpy_arrays_in_csv_empty_dataframe()
    test_preserve_numpy_arrays_in_csv_mixed_dtypes()
    print("\nDone.")
