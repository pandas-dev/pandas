import math

a = 1450804465901089690
b = 1450804465901089614

print("difference:", a - b)
print("float(a):", float(a))
print("float(b):", float(b))
print("float difference:", float(a) - float(b))
print(math.isclose(float(a), float(b), rel_tol=0, abs_tol=100))