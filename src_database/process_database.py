# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    process_database.py                                :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/11 17:43:17 by Danilo            #+#    #+#              #
#    Updated: 2023/07/11 17:51:04 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

#    1) Starts the event 5 days before the notice date and ends it 2 days later                                                           
#    2) Removes duplicate dates from the original database
#    3) If there are events at least 7 days apart, combines them into one event
#    4) Filters out rows with dates older than 2020 (Waverys limit)

import pandas as pd 

df = pd.read_csv('../database/noticias-ressacas_v2.csv', encoding='latin1')
dates = pd.to_datetime(df['Data'], dayfirst=True)

# Calculate start and end dates
start_dates = dates - pd.DateOffset(days=5)
end_dates = dates + pd.DateOffset(days=2)

# Create new DataFrame with 'start' and 'end' columns
df = pd.DataFrame({'dates': dates, 'start': start_dates, 'end': end_dates})

# Remove duplicate data
df = df.drop_duplicates()

# Group the DataFrame by dates
grouped = df.groupby('dates')

# Define a function to check for dates at least 7 days apart and combine start and end
def combine_dates(group):
    if len(group) > 1 and (group['dates'].max() - group['dates'].min()).days >= 7:
        start = group['start'].min()
        end = group['end'].max()
        return pd.Series({'start': start, 'end': end})
    else:
        return group.iloc[0][['start', 'end']]

# Apply the function to each group and create a new DataFrame
new_df = grouped.apply(combine_dates).reset_index(drop=True)

# Convert 'start' column to the desired format (yyyymmdd)
new_df['id'] = new_df['start'].dt.strftime('%Y%m%d')

# Filter out rows with dates older than 2020
new_df = new_df[new_df['start'].dt.year >= 2020]

# Reset the index of the DataFrame
new_df = new_df.reset_index(drop=True)

# Set 'id' as the index
new_df = new_df.set_index('id')

# Save the DataFrame
new_df.to_csv('../database/dates-limits.csv')