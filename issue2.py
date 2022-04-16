import pandas as pd
import json

df = pd.DataFrame({
    "X": [1,2],
    "Y": [1,2]}
)
for x, df in df.groupby("X"):
    print(type(x))
    print(x)  # <class 'int'>

for (x, y), df in df.groupby(["X","Y"]):
    print(int(x))
    print(type(x))
    # print(y)

    update_data = {
         'x_value': int(x),
         'y_value': int(y)
         }
    stringtype = json.dumps(update_data)
    print(stringtype)

    


#for ur understanding
# ipl_data = {'Team': ['Riders', 'Riders', 'Devils', 'Devils', 'Kings',
#    'kings', 'Kings', 'Kings', 'Riders', 'Royals', 'Royals', 'Riders'],
#    }
# df = pd.DataFrame(ipl_data)

# print (df.groupby('Team').groups)






   
