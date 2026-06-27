# clean.py
''' 
This script will import and clean the positional data from four ground stations along the 
cascades region.
Jun 20, 2026
Version 1
Chelsea Momoh
'''

import pandas as pd
import json

#Step 1: Document the Coordinate Reference Frame
metadata = {}

#Filepaths is a dictionary with key:value pairs "filename":"filepath"
filepaths = {
    "P349": "/workspaces/GNSS/data/P349.cwu.nam14.csv",
    "P380": "/workspaces/GNSS/data/P380.cwu.nam14.csv",
    "P434": "/workspaces/GNSS/data/P434.cwu.nam14.csv",
    "P441": "/workspaces/GNSS/data/P441.cwu.nam14.csv"
}

print("Beginning Data Clean.")
#For key, value in the filepaths dictionary
for station_id, path in filepaths.items():
    #Open a csv file
    with open(path, 'r') as f:
        #Take the first ten lines, strip them of blankspace, and store them in a list object.
        lines = [next(f).strip() for _ in range(10)]
    #Now add key vaue pairs to the metadata dict for later reference 
    metadata[station_id] = lines
    #To confirmt he type of object that lines is
    #print(type(lines))
    #print(type(lines[0]))

#I'd liek to save the data persistently
with open("metadata.json", "w") as f:
    json.dump(metadata, f, indent=4)


#Step 2: Handle missing dsata

'''print("\n\nCatching missing values.")
missing = p349.isnull().sum()
print(missing)
#print((p349 == 'NaN').sum()) This is not a missing calues format

missing = p380.isnull().sum()
print(missing)

missing = p434.isnull().sum()
print(missing)

missing = p441.isnull().sum()
print(missing)'''

print("\nAfter previous data handling, I found there were no missing values in the dataset\n")


#Step 3: Change date column type to datetime
p349 = pd.read_csv("/workspaces/GNSS/data/P349.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p380 = pd.read_csv("/workspaces/GNSS/data/P380.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p434 = pd.read_csv("/workspaces/GNSS/data/P434.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')
p441 = pd.read_csv("/workspaces/GNSS/data/P441.cwu.nam14.csv", skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')

'''
#CAREFUL! THIS CREATES COPIES AND DOES NOT MODIFY THE ORIGINAL. Pandas!
datasets = [("P349", p349), ("P380", p380), ("P434", p434), ("P441", p441)]
for station_id, dataset in datasets:
    dataset['Date'] = pd.to_datetime(dataset['Date'])
'''

p349['Date'] = pd.to_datetime(p349['Date'])
p380['Date'] = pd.to_datetime(p380['Date'])
p434['Date'] = pd.to_datetime(p434['Date'])
p441['Date'] = pd.to_datetime(p441['Date'])

#Step 4: creating JSON files for persistent storage 
p349.to_json("p349.json", orient="records", date_format="iso", indent=4)
p380.to_json("p380.json", orient="records", date_format="iso", indent=4)
p434.to_json("p434.json", orient="records", date_format="iso", indent=4)
p441.to_json("p441.json", orient="records", date_format="iso", indent=4)
print("Completed Data Cleaning")


#Never running this code again. 15,000 lines takes a toll on my computer apparently.
