import pandas as pd

values = [
    {"col1": 30.0, "col2": 116.80000305175781},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
    {"col1": 30.100000381469727, "col2": 116.8000030517578},
    {"col1": None, "col2": None},
    {"col1": None, "col2": None},
]

data = pd.DataFrame(values)
print(data.corr(method="pearson"))

# def corr_coef(X,Y):
#     x1 = np.array(X)
#     y1 = np.array(Y)
#     x_m=x1.mean()
#     y_m=y1.mean()
#     number=0
#     v1sq=0
#     v2sq=0
#     for i in range(len(x1)):
#         xx = (x1[i]-x_m)
#         yy = (y1[i]-y_m)
#         number+=xx*yy
#         v1sq+=xx*xx
#         v2sq+=yy*yy
#     return(number/(math.sqrt(v1sq*v2sq)))

data = data.dropna()
# print(corr_coef(data.iloc[:,0],data.iloc[:,1]))
