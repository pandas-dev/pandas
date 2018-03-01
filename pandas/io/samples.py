import numpy
import pandas


def _countries_with_penguins():
    """Main dataset to be used in pandas documentation examples."""
    columns = ('Code', 'Name', 'Capital', 'Continent',
               'Penguin species', 'Avg. temperature')
    data = [
        ('AO', 'Angola', 'Luanda', 'AF', 1, 21.55),
        ('AQ', 'Antarctica', numpy.nan, 'AN', 7, numpy.nan),
        ('AR', 'Argentina', 'Buenos Aires', 'SA', 7, 14.8),
        ('AU', 'Australia', 'Canberra', 'OC', 10, 21.65),
        ('BR', 'Brazil', 'Brasilia', 'SA', 1, 24.95),
        ('CL', 'Chile', 'Santiago', 'SA', 10, 8.45),
        ('EC', 'Ecuador', 'Quito', 'SA', 2, 21.85),
        ('MZ', 'Mozambique', 'Maputo', 'AF', 1, 23.8),
        ('NA', 'Namibia', 'Windhoek', 'AF', 2, 19.95),
        ('NZ', 'New Zealand', 'Wellington', 'OC', 7, 10.55),
        ('PE', 'Peru', 'Lima', 'SA', 2, 19.6),
        ('UY', 'Uruguay', 'Montevideo', 'SA', 3, 17.55),
        ('ZA', 'South Africa', 'Pretoria', 'AF', 4, 17.75),
    ]
    df = pandas.DataFrame(data, columns=columns).set_index('Code')
    df['Continent'] = df['Continent'].astype('category')
    return df


def _stocks():
    """Sample dataset to be used in time series examples in pandas docs."""
    data = {'MSFT': [62.9423, 62.2900, 64.5091, 67.0557, 68.4074, 67.9031,
                     71.6169, 73.6561, 73.7712, 82.3774, 83.3578, 85.1401],
            'RHT': [75.8800, 82.8100, 86.5000, 88.0800, 89.5700, 95.7500,
                    98.8700, 107.5000, 110.8600, 120.8300, 126.7600, 120.1000],
            'AAPL': [118.9445, 134.2745, 141.4228, 141.4129, 150.3810,
                     142.3624, 147.0181, 162.1124, 152.9444, 167.7505,
                     170.5391, 168.5428]}
    index = pandas.to_datetime(['2017-01-01', '2017-02-01', '2017-03-01',
                                '2017-04-01', '2017-05-01', '2017-06-01',
                                '2017-07-01', '2017-08-01', '2017-09-01',
                                '2017-10-01', '2017-11-01', '2017-12-01'])
    index.name = 'Date'
    df = pandas.DataFrame(data, index=index)
    return df


countries_with_penguins = _countries_with_penguins()
stocks = _stocks()
