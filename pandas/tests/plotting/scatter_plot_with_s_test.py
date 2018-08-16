import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

#n = 50

#data = np.random.rand(n, 3)
#df = pd.DataFrame(data)
#df.columns = ['a', 'b', 'c']
#sizes = np.random.rand(n) * 100.0

#df.plot.scatter(x='a', y='b', s='c')

hp = [100, 104.4, 85, 180, 220, 100, 115, 500, 350, 410,
      530, 550, 480, 100, 120, 105, 90, 108, 130, 80]

ms = [180, 203, 170, 255, 280, 187, 195, 170, 145, 160,
      170, 180, 165, 280, 300, 270, 245, 290, 310, 240]

vehicle_type = ['Car', 'Car', 'Car', 'Car', 'Car', 'Car', 'Car',
                'Truck', 'Truck', 'Truck', 'Truck', 'Truck', 'Truck',
                'Motorcycle', 'Motorcycle', 'Motorcycle', 'Motorcycle',
                'Motorcycle', 'Motorcycle', 'Motorcycle']

cat_type= CategoricalDtype(categories=['Motorcycle', 'Car', 'Truck'],
                           ordered=True)
df = pd.DataFrame()
df['Horse power'] = hp
df['Max speed'] = ms
df['Vehicle type'] = pd.Series(vehicle_type).astype(cat_type)
ax = df.plot.scatter(x='Horse power', y='Max speed', s='Vehicle type')

plt.show()
