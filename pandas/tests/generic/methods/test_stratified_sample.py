import pandas as pd
import pandas._testing as tm


def test_stratified_sample():
    df = pd.DataFrame({
        'gender': ['male', 'male', 'female', 'female', 'female',
                   'female', 'female', 'female', 'male', 'male'],
        'age': [25, 26, 25, 26, 30, 25, 25, 30, 30, 25],
        'country': ['US', 'CAN', 'MEX', 'CAN', 'IN', 'CAN', 'CAN', 'US', 'CAN', 'IN'],
        'income_K' : [100, 110, 99, 110, 110, 100, 100, 110, 100, 99]
    })

    result = df.stratified_sample(strata=['age'], n=5, random_state=0, reset_index=True)
    expected = pd.DataFrame({
        'gender': ['female', 'male', 'female', 'male', 'female'],
        'age': [25, 25, 26, 30, 30],
        'country': ['CAN', 'US', 'CAN', 'CAN', 'US'],
        'income_K' : [100, 100, 110, 100, 110]
    })

    tm.assert_frame_equal(result, expected)


def test_stratified_sample_counts():
    df = pd.DataFrame({
        'gender': ['male', 'male', 'female', 'female', 'female',
                   'female', 'female', 'female', 'male', 'male'],
        'age': [25, 26, 25, 26, 30, 25, 25, 30, 30, 25],
        'country': ['US', 'CAN', 'MEX', 'CAN', 'IN', 'CAN', 'CAN', 'US', 'CAN', 'IN'],
        'income_K' : [100, 110, 99, 110, 110, 100, 100, 110, 100, 99]
    })

    result = df.stratified_sample_counts(strata=['age'], n=5)
    expected = pd.DataFrame({
        'age': [25, 26, 30],
        'size': [5, 2, 3],
        'sample_size': [2, 1, 2]
    })

    tm.assert_frame_equal(result, expected)
