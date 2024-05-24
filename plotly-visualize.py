import pandas as pd
import plotly.express as px

def plotly_visualize(df, x, y, chart_type='scatter', title=None):
    """
    Create an interactive plot using Plotly from a Pandas DataFrame.
    
    Parameters:
    - df: Pandas DataFrame
    - x: column name for x-axis
    - y: column name for y-axis
    - chart_type: type of the chart ('scatter', 'line', 'bar')
    - title: title of the chart
    """
    if chart_type == 'scatter':
        fig = px.scatter(df, x=x, y=y, title=title)
    elif chart_type == 'line':
        fig = px.line(df, x=x, y=y, title=title)
    elif chart_type == 'bar':
        fig = px.bar(df, x=x, y=y, title=title)
    else:
        raise ValueError(f"Unsupported chart type: {chart_type}")
    
    fig.show()

data = {
    'A': [1, 2, 3, 4, 5],
    'B': [5, 4, 3, 2, 1]
}
df = pd.DataFrame(data)
# plotly_visualize(df, x='A', y='B', chart_type='scatter', title='Scatter Plot Example')