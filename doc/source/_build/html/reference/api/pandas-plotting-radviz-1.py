df = pd.DataFrame(
    {
        'SepalLength': [6.5, 7.7, 5.1, 5.8, 7.6, 5.0, 5.4, 4.6, 6.7, 4.6],
        'SepalWidth': [3.0, 3.8, 3.8, 2.7, 3.0, 2.3, 3.0, 3.2, 3.3, 3.6],
        'PetalLength': [5.5, 6.7, 1.9, 5.1, 6.6, 3.3, 4.5, 1.4, 5.7, 1.0],
        'PetalWidth': [1.8, 2.2, 0.4, 1.9, 2.1, 1.0, 1.5, 0.2, 2.1, 0.2],
        'Category': [
            'virginica',
            'virginica',
            'setosa',
            'virginica',
            'virginica',
            'versicolor',
            'versicolor',
            'setosa',
            'virginica',
            'setosa'
        ]
    }
)
pd.plotting.radviz(df, 'Category')
