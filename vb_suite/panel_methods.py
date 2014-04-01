from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# shift

setup = common_setup + """
index = date_range(start="2000", freq="D", periods=1000)
panel = Panel(np.random.randn(100, len(index), 1000))
"""

panel_shift = Benchmark('panel.shift(1)', setup,
                               start_date=datetime(2012, 1, 12))

panel_shift_minor = Benchmark('panel.shift(1, axis=minor)', setup,
                               start_date=datetime(2012, 1, 12))
