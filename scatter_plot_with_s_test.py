import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#n = 50
#data = np.random.rand(n, 4)
#df = pd.DataFrame(data)
#df.columns = ['a', 'b', 'c', 'd']
#sizes = np.random.rand(n) * 100.0
#df.plot.scatter(x='a', y='b', s='c', c='d', colormap='RdYlGn')

surf_area = [(2*x+1.0)*40 for x in np.random.rand(30)]+ \
            [(2*x+1.0)*80 for x in np.random.rand(20)]+ \
            [(2*x+1.0)*100 for x in np.random.rand(10)]

surf_area = np.array(surf_area)

types = ['Flat' for i in range(30)] + ['House' for i in range(20)] + ['Castle' for i in range(10)]

types = np.array(types)

prices = 0.01 * surf_area * (np.random.rand(60) + 1.5) / 2
prices[30:50] *= 1.4
prices[50:] *= 2

df = pd.DataFrame()
df['Surface area (sqm)'] = surf_area
df['Price (M€)'] = prices
df['Property type'] = pd.Categorical(types, categories=['Flat', 'House', 'Castle'], ordered=True)

df.plot.scatter(x='Surface area (sqm)', y='Price (M€)',
                s='Property type', alpha=.5)

#hp = [100, 104.4, 85, 180, 220, 100, 115, 500, 350, 410,
#      530, 550, 480, 100, 120, 105, 90, 108, 130, 80]

#ms = [180, 203, 170, 255, 280, 187, 195, 170, 145, 160,
#      170, 180, 165, 280, 300, 270, 245, 290, 310, 240]

#df = pd.DataFrame()
#df['Horse power'] = hp
#df['Max speed'] = ms

#vehicle_type = ['Car', 'Car', 'Car', 'Car', 'Car', 'Car', 'Car',
#                'Truck', 'Truck', 'Truck', 'Truck', 'Truck', 'Truck',
#                'Motorcycle', 'Motorcycle', 'Motorcycle', 'Motorcycle',
#                'Motorcycle', 'Motorcycle', 'Motorcycle']
#categories=['Motorcycle', 'Car', 'Truck']
#df['Vehicle type'] = pd.Categorical(vehicle_type, categories=categories, ordered=True)
#ax = df.plot.scatter(x='Horse power', y='Max speed', s='Vehicle type', size_factor=.5)
plt.show()