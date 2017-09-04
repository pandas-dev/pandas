import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.DataFrame(np.random.rand(20, 4), columns = ['Width' ,'Height', 'Weight', 'Price'])

df['Price'] = np.arange(1, 21)

df1 = pd.DataFrame({
    'Width': np.arange(10),
    'Height': np.ones(10),
    'Price': [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    })

fig1, ax1 = plt.subplots(figsize = (5, 5))
fig2, ax2 = plt.subplots(figsize = (5, 5))

df.plot(kind="scatter", x='Width', y='Height', s='Price', s_grow= 10, ax=ax1)
df.plot(kind="scatter", x='Width', y='Height', s='Price', s_grow= 1, ax=ax2)

plt.show()

