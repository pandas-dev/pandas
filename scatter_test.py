import matplotlib.pyplot as plt
import pandas as pd

midwest = pd.read_csv("http://goo.gl/G1K41K")
midwest = midwest[midwest['poptotal'] < 50000]
graph = midwest.plot(kind='scatter',
                     x='area',
                     y='poptotal',
                     s='popdensity',
                     s_grow=0.2,
                     title='Popuation vs area and density')
plt.show()
