import pandas as pd

# Acceptance test for issue 43629
print("Create and validate an interval")
interval = pd.Interval(1,2,"neither")
print("Interval left hand side:")
print(interval.left)
print("Expected Result:")
print(1)
print("Interval right hand side:")
print(interval.right)
print("Expected Result:")
print(2)
print("Interval closed sides:")
print(interval.closed)
print("Expected Result:")
print("neither")
print('++++++++++++++++++++++++++++++++++++++++++++')

print("Add two intervals and validate their sum")
interval_a = pd.Interval(0, 10)
interval_b = pd.Interval(1, 2)
interval_c = interval_a + interval_b
print("Actual Result:")
print(interval_c)
print("Expected Result:")
print(pd.Interval(0, 10))
print('++++++++++++++++++++++++++++++++++++++++++++')

print("Subtract two intervals and validate their difference")
interval_d = interval_c - interval_a
print("Actual Result:")
print(interval_d)
print("Expected Result:")
print(pd.Interval(0, 0,'neither'))
print('++++++++++++++++++++++++++++++++++++++++++++')

print("Multiply two intervals and validate their product")
interval_e = interval_c * interval_b
print("Actual Result:")
print(interval_e)
print("Expected Result:")
print(pd.Interval(1, 2))
print('++++++++++++++++++++++++++++++++++++++++++++')

# Create an interval array of the previous created intervals a and b
print("Create and validate an interval array")
interval_array = pd.arrays.IntervalArray([interval_a, interval_b])
print("Actual Result:")
print(interval_array)
print("Expected Result:")
print(
'''<IntervalArray>
[(0, 10], (1, 2]]
Length: 2, dtype: interval[int64, right]'''
)

print("Calculate and verify the sum of the interval array")
interval_arr_sum = interval_array.sum()
print("Actual Result:")
print(interval_arr_sum)
print("Expected Result:")
print(pd.arrays.IntervalArray([pd.Interval(0, 10)]))
print('++++++++++++++++++++++++++++++++++++++++++++')

print("Calculate and verify the product of the interval array")
interval_arr_product = interval_array.product()
print("Actual Result:")
print(interval_arr_product)
print(pd.arrays.IntervalArray([pd.Interval(1, 2)]))
print('++++++++++++++++++++++++++++++++++++++++++++')