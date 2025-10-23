import pandas as pd
df = pd.DataFrame({'a' : [1, 2, 3, 4], 'b' : [4, 3, 2, 1]})
df['b'] = df['b'].astype('category').cat.set_categories([4, 3, 2, 1], ordered=True)
#import pdb; pdb.set_trace()
crr = df.corr(method='spearman')
print(crr)


df = pd.DataFrame({'a' : [1, 2, 3, 4], 'b' : ["vh", "h", "m", "l"]})
df['b'] = df['b'].astype('category').cat.set_categories(["vh", "h", "m", "l"], ordered=True)
#import pdb; pdb.set_trace()
print(df)
print(df.dtypes)
crr = df.corr(method='spearman')
print(crr)

ser_ord_cat = pd.Series( pd.Categorical(
             ["vh", "h", "m", "low"], 
             categories=["vh", "h", "m", "low"], ordered=True
             ))
print(ser_ord_cat)
crr = ser_ord_cat.corr(ser_ord_cat, method='spearman')
print(crr)