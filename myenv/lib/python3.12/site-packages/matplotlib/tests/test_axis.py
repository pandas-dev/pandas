import numpy as np

import matplotlib.pyplot as plt
from matplotlib.axis import XTick


def test_tick_labelcolor_array():
    # Smoke test that we can instantiate a Tick with labelcolor as array.
    ax = plt.axes()
    XTick(ax, 0, labelcolor=np.array([1, 0, 0, 1]))
