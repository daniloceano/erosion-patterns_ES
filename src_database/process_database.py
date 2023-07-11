import pandas as pd 

df = pd.read_csv('../database/noticias-ressacas_v2.csv', encoding='latin1')
dates = pd.to_datetime(df['Data'], dayfirst=True)
pd.Series(dates.unique())