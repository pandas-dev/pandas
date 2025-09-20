from numba import typed, int64

# used in TestMiscErrorHandling::test_handling_of_write_to_*_global
_global_list = [1, 2, 3, 4]

_global_dict = typed.Dict.empty(int64, int64)


def global_reflected_write():
    _global_list[0] = 10


def global_dict_write():
    _global_dict[0] = 10
