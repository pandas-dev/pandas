import pandas as pd;
print("version",pd.__version__)

# date_interval = pd.interval_range(pd.Timestamp('2018-01-01'), freq='3D', periods=3)
# print("data",date_interval)
# table_format = pd.DataFrame({'date': date_interval})
# print("Table Format")
# print(table_format)

# pi = pd.period_range('2018-01-01', freq='D', periods=3)
# print("sequence of date", pi)

# print("different between list and astype")
# print("list of object--> ",list(pi), " only date--> ", list(pi.astype(str)))
# array1 = date_interval.get_indexer(list(pi.astype(str)))
# print(array1)

date_interval = pd.interval_range(pd.Timestamp('2018-01-01'), freq='3D', periods=3)
print("interval range", date_interval)
pi = pd.period_range('2018-01-05', freq='D', periods=3)
print("value to test", pi)
array = date_interval.get_indexer(pi)
print("array of object",pd.array(pi))
array1 = date_interval.get_indexer(pd.array(pi))
print("list of object",list(pi))
array3 = date_interval.get_indexer(list(pi))
print("value to test",list(pi) + [10])
array4 = date_interval.get_indexer(list(pi) + [10])
print("wrong output is", array4)

print("correct format",list(pi.astype(str)))
output = date_interval.get_indexer(list(pi.astype(str)))
print("Proper output", output)




