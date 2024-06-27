from numba import int32
from numba.typed import List


# global typed-list for testing purposes
global_typed_list = List.empty_list(int32)
for i in (1, 2, 3):
    global_typed_list.append(int32(i))


def catch_global():
    x = List()
    for i in global_typed_list:
        x.append(i)
