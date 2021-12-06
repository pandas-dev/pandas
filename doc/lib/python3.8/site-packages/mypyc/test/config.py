import os

this_file_dir = os.path.dirname(os.path.realpath(__file__))
prefix = os.path.dirname(os.path.dirname(this_file_dir))

# Locations of test data files such as test case descriptions (.test).
test_data_prefix = os.path.join(prefix, 'mypyc', 'test-data')
