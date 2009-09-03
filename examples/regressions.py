from datetime import datetime
import string

from pandas.core.api import DataMatrix, DateRange

N = 100

start = datetime(2009, 9, 2)
dateRange = DateRange(start, periods=N)

def makeDataMatrix():
    data = DataMatrix(np.random.randn(N, 7),
                      columns=list(string.ascii_uppercase[:7]),
                      index=dateRange)

    return data

data = makeDataMatrix()
model = LinearModel(data, window=50, minPeriods=20)
model.parseFormula('A ~ B + C + D + E + F + G + I')
model.fit()

data = {
    'A' : makeDataMatrix(),
    'B' : makeDataMatrix(),
    'C' : makeDataMatrix()
}

panelModel = XSLinearModel(data, window=50, minPeriods=20)
panelModel.parseFormula('A ~ B + C + I')
panelModel.fit()
