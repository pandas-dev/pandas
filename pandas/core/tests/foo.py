from test_sparse import *

data_dict = {
    'item1' : panel_data1(),
    'item2' : panel_data2(),
    'item3' : panel_data3()
}
panel = SparseWidePanel(data_dict)

