# clean.py
''' 
This script will import and clean the positional data from four ground stations along the 
cascades region.
Jun 20, 2026
Version 1
Chelsea Momoh
'''

import pandas as pd

#Step 0: Import Data and as Dataframe
#Each row in th ecsv has trailing comma, so I use lambda to delete the last row, which was initially called "Unnamed: 8"
p349 = pd.read_csv("/workspaces/GNSS/data/P349.cwu.nam14.csv", 
    skipinitialspace=True,
    usecols=lambda x: x != 'Unnamed: 8')
p380 = pd.read_csv("/workspaces/GNSS/data/P380.cwu.nam14.csv",
    skipinitialspace=True,
    usecols=lambda x: x != 'Unnamed: 8')
p434 = pd.read_csv("/workspaces/GNSS/data/P434.cwu.nam14.csv", 
    skipinitialspace=True,
    usecols=lambda x: x != 'Unnamed: 8')
p441 = pd.read_csv("/workspaces/GNSS/data/P441.cwu.nam14.csv", 
    skipinitialspace=True,
    usecols=lambda x: x != 'Unnamed: 8')




#Step 1: Document the Coordinate Reference Frame
metadata = {}

datasets = [("P349", p349), ("P380", p380), ("P434", p434), ("P441", p441)]

for station_id, dataset in datasets:
    metadata[station_id] = {
        "Format Version": dataset.iloc[1].to_list(),
        "Reference Frame": dataset.iloc[2].to_list(),
        "4-character ID": dataset.iloc[3].to_list(),
        "Station name": dataset.iloc[4].to_list(),
        "Begin Date": dataset.iloc[5].to_list(),
        "End Date": dataset.iloc[6].to_list(),
        "Release Date": dataset.iloc[7].to_list(),
        "Source File": dataset.iloc[8].to_list(),
        "Offset from source file": dataset.iloc[9].to_list(),
        "Reference position": dataset.iloc[10].to_list()
    }


print("Printing keys for each dataset's metadata dictionary.")    
for key in metadata.keys():
    print(key)

print("\n\nPrinting P349 metadata:")
print(metadata['P349'])

# Deleting metadata from dataframe
p349 = pd.read_csv("/workspaces/GNSS/data/P349.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p380 = pd.read_csv("/workspaces/GNSS/data/P380.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p434 = pd.read_csv("/workspaces/GNSS/data/P434.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p441 = pd.read_csv("/workspaces/GNSS/data/P441.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')

datasets = [("P349", p349), ("P380", p380), ("P434", p434), ("P441", p441)]
for station_id, dataset in datasets:
    dataset['Date'] = pd.to_datetime(dataset['Date'])

'''print("\n\nPrinting dataset head")
print(p349.head())
print(p349.columns.tolist())'''


#Step 2: Handle missing dsata

print("\n\nCatching missing values.")
missing = p349.isnull().sum()
print(missing)
#print((p349 == 'NaN').sum()) This is not a missing calues format

missing = p380.isnull().sum()
print(missing)

missing = p434.isnull().sum()
print(missing)

missing = p441.isnull().sum()
print(missing)

print("There are no missing values in the dataset")









'''with open("/workspaces/GNSS/data/P349.cwu.nam14.csv", 'r') as f:
    for i, line in enumerate(f):
        if i < 11:
            metadata[i] = line.strip()
        else:
            break'''


            #Next: must properly load datasets to correct Iloc issue