
date_index = pd.date_range('1/1/2010', periods=6, freq='D')
df2 = pd.DataFrame({"prices": [100, 101, np.nan, 100, 89, 88]},
                   index=date_index)
date_index2 = pd.date_range('12/29/2009', periods=10, freq='D')
df2=df2.reindex(date_index2)
