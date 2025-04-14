import pandas as pd

pd.set_option('display.float_format', '{:.2f}'.format)

df = pd.DataFrame([[15.22345676543234567]*6], columns=[1,2,3,4,5,6])

# Default float_format works:
print(df)  # ✅ Columns display with 2 decimals

default_fmt = '{:.2f}'
special_fmt = {1: '{:.1f}'}

formats = {
    col: special_fmt.get(col, default_fmt)
    for col in df.columns
}

styled = df.style.format(formats)
styled.to_html("styled.html")
