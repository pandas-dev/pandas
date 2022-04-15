import pandas as pd


def test_group_last_use_ops_0():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype() 


def test_group_last_use_ops_1():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_2():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_3():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_4():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_5():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_6():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_6():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_7():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_8():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()

def test_group_last_use_ops_9():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_10():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_11():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_12():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_13():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_14():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_15():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_16():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_17():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_18():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_19():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, pd.NA, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_20():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, pd.NA, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_21():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, pd.NA, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_22():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, pd.NA, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_23():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, pd.NA, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_24():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, pd.NA, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()



def test_group_last_use_ops_26():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_27():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_28():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, True, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_29():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, False, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_30():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, True, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_31():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, False, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_32():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, pd.NA, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_33():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, pd.NA, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_34():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, True, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_35():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, False, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_36():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, pd.NA, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_37():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, pd.NA, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_38():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, pd.NA, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_39():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, pd.NA, pd.NA, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_40():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_41():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_42():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, True, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_43():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, False, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_44():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, True, True, pd.NA, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_45():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, False, False, pd.NA, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_46():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, True, True, True, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_47():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, False, False, False, pd.NA]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_48():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, pd.NA, pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_49():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, pd.NA, pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_50():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, pd.NA, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_51():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, pd.NA, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_52():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, pd.NA, pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_53():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, pd.NA, pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_55():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [True, True, pd.NA, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_56():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [False, False, pd.NA, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_57():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [True, pd.NA, pd.NA, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_58():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [False, pd.NA, pd.NA, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_59():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [True, True, pd.NA, pd.NA, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_60():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [False, False, pd.NA, pd.NA, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_61():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [pd.NA, pd.NA, pd.NA, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_62():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [pd.NA, pd.NA, pd.NA, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_63():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_64():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_65():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_66():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, True, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_67():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, False, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_68():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, True, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_69():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, False, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_70():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, False, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_71():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_72():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, False, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_73():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, True, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_74():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, True, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_75():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_76():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, False, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_77():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, False, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_78():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_79():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_80():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_81():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [True, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_82():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [False, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_83():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, pd.NA, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_84():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, pd.NA, True, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_85():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, pd.NA, False, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_86():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4],
            'test': [pd.NA, pd.NA, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_87():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_88():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_89():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_90():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_91():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [True, pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_92():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [False, pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_93():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, pd.NA, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_94():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5],
            'test': [pd.NA, pd.NA, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_95():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_96():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_97():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_98():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_99():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [True, pd.NA, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_100():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [False, pd.NA, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()

def test_group_last_use_ops_101():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, pd.NA, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_102():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6],
            'test': [pd.NA, pd.NA, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_103():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [True, True, True, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_104():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [False, False, False, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_105():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [pd.NA, True, True, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_106():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [pd.NA, False, False, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_107():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [True, pd.NA, True, True, True, True, True, True, True, True, True, True,
                     True, True, True, True, True, True, True, True]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_108():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [False, pd.NA, False, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()


def test_group_last_use_ops_109():
    df = pd.DataFrame(
        {
            'id': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            'test': [pd.NA, pd.NA, False, False, False, False, False, False, False, False,
                     False, False, False, False, False, False, False, False, False, False]
        }
    ).convert_dtypes()
    grouped = df.groupby('id')
    bad = grouped.last()
    
    assert bad.test.dtype == pd.BooleanDtype()

