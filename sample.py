import pandas as pd
import numpy as np
import csv

df = pd.DataFrame({"col": [8.5557]}, dtype=np.float32)

print(df.to_csv())

print(df.to_csv(quoting=csv.QUOTE_NONNUMERIC))

values = np.array([8.57], dtype="float32")
print(values)
# [8.57]

print(np.array(values, dtype="object"))
# [8.569999694824219]

print(np.array(values, dtype="str"))


# Original array in float32
float32_arr = np.array([1.2345678, 2.3456789], dtype=np.float32)

# Convert to object
object_arr = float32_arr.astype(object)

print("Original float32 array:", float32_arr)
print("Object array:", object_arr)
print("Data type of object_arr:", object_arr.dtype)
