# script to generate the pandas logo

from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np

rcParams['mathtext.fontset'] = 'cm'


def fnx():
    return np.random.randint(5, 50, 10)


fig = plt.figure(figsize=(6, 1.25))

ax = fig.add_axes((0.45, 0.1, 0.16, 0.8))
bar_data = [2.1, -00.8, 1.1, 2.5, -2.1, -0.5, -2.0, 1.5]
ax.set_ylim(-3, 3)
ax.set_xticks([])
ax.set_yticks([])
ax.bar(np.arange(len(bar_data)), bar_data)

ax = fig.add_axes((0.63, 0.1, 0.16, 0.8))
for i in range(4):
    ax.plot(np.random.rand(8))
ax.set_xticks([])
ax.set_yticks([])

ax = fig.add_axes((0.63 + 0.18, 0.1, 0.16, 0.8))
y = np.row_stack((fnx(), fnx(), fnx()))
x = np.arange(10)
y1, y2, y3 = fnx(), fnx(), fnx()
ax.stackplot(x, y1, y2, y3)
ax.set_xticks([])
ax.set_yticks([])

plt.figtext(0.05, 0.5, "pandas", size=40)

plt.figtext(
    0.05, 0.2, r"$y_{it} = \beta^{\prime} x_{it} + \mu_{i} + \epsilon_{it}$",
    size=16, color="#5a89a4")

fig.savefig('pandas_logo.svg')
fig.savefig('pandas_logo.png')
