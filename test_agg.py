import pandas as pd
import numpy as np

def f(x):
   return list(x)

#df = pd.DataFrame({'A' : [1, 1, 3], 'B' :  [1, 2, 4]})
#result = df.groupby('A').aggregate(f)


#df = pd.DataFrame({'A' : [1, 1, 3], 'B' :  [1, 2, 4]})
#result = df.groupby('A').aggregate(list)
#result = df.groupby('A').agg(list)

df = pd.DataFrame({'A' : [1, 1, 3], 'B' :  [1, 1, 4], 'C' :  [1, 3, 4]})
#result = df.groupby(['A', 'B']).aggregate(pd.Series)


#df = pd.DataFrame({'A': [1, 1, 1, 3, 3, 3],
 #                          'B': [1, 1, 1, 4, 4, 4], 'C': [1, 1, 1, 3, 4, 4]})

#print ('series ')
result = df.groupby('A')['C'].aggregate(np.array)
#print (result)
#
result = df.groupby(['A', 'B']).aggregate(np.array)
#print (result)
#
# result = df.groupby('A')['C'].aggregate(list)
# print (result)

def f(x):
   return np.array(x)

print ('array')
result = df.groupby(['A', 'B']).aggregate(f)
print (result)

# result = df.groupby('A')['C'].aggregate(tuple)
# expected = pd.Series([(1, 1, 1), (3, 4, 4)], index=[1, 3], name='C')
# expected.index.name = 'A'



