import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

x = [10,50,40]
y = "male","female","unknown"

data = [['male', 10], ['female', 50], ['unknown', 40]] 
  
df = pd.DataFrame(data, columns = ['x', 'y']) 


plt = df.plot.scatter(x='x', y='y')
plt.set(xlabel="gender", ylabel = "age")