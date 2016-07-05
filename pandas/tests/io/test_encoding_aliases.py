import pandas, os, nose

def test_read_csv():
    # see gh issue 13549
    test_encodings = ['utf-16', 'utf_16', 'UTF_16', 'UTF-16']
    engines = ['c', 'python', None]
    expected_output = [1]*len(test_encodings)*len(engines)
    test_output = []
    path = "test.csv"
    pandas.DataFrame({"A": [0,1], "B": [2,3]}).to_csv(
        path, encoding="utf-16")

    for encoding in test_encodings:
        for engine in engines:
            try:
                pandas.io.parsers.read_csv(
                    path,
                    engine='c',
                    encoding=encoding)
                print(encoding, 'succeeded with engine =', engine)
                test_output.append(1)
            except UnicodeDecodeError:
                print(encoding, 'failed with engine =', engine)
                test_output.append(0)

    assert (expected_output == test_output)

    os.remove("test.csv")

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)