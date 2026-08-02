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

# Document the Coordinate Reference Frame
metadata = {}

#Filepaths is a dictionary with key:value pairs "filename":"filepath"
filepaths = {
    "P349": "/workspaces/GNSS/data/P349.cwu.nam14.csv",
    "P380": "/workspaces/GNSS/data/P380.cwu.nam14.csv",
    "P434": "/workspaces/GNSS/data/P434.cwu.nam14.csv",
    "P441": "/workspaces/GNSS/data/P441.cwu.nam14.csv"
}


class Station:
    def __init__(self, station_id, filepath):
        #These are attributes of the station class. They will become attributes of each station object.
        self.station_id = station_id 
        self.filepath = filepath
        self.metadata = self.extract_metadata()
        self.dataframe = pd.read_csv(filepath, skiprows=11, skipinitialspace=True, usecols=lambda x: x != 'Unnamed: 8')

    def extract_metadata(self):
        #These are methods of the station class. Functions that each station object can perform.
        with open(self.filepath, 'r') as f:
            lines = [next(f).strip() for _ in range(10)]
        return lines

    def handle_missing_values(self):
        # There were no missing values on my original run, so I didn't do anything to handle them. 
        # It would be better, in the future, to implement a handling method.
        missing = self.dataframe.isnull().sum()
        print(f"Missing values for {self.station_id}:")
        print(missing)
        return missing 

    def change_date_column_type(self):
        self.dataframe['Date'] = pd.to_datetime(self.dataframe['Date'])
        return self.dataframe

    def save_to_json(self):
        json_filename = f"{self.station_id}.json"
        self.dataframe.to_json(json_filename, orient="records", date_format="iso", indent=4)
        print(f"Saved {self.station_id} data to {json_filename}")

    def clean_data(self):
        self.handle_missing_values() #These methods don't require arguments. They just utilize attributes of the class object.
        self.dataframe = self.change_date_column_type()
        self.save_to_json()

print("Beginning Data Clean.")

for station_id, path in filepaths.items():
    station = Station(station_id, path)
    metadata[station_id] = station.metadata
    station.clean_data()





    # ------- OUTDATED //// PRIOR to REFACTORING //// OUTDATED -------

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

# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
#Note on project limitation: Sampling Gaps  
# * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -
# I assumed that a measurement was taken each day without checking for "skipped" days. 
# This does not negatively affect my results, but it is better to not skip this step.


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
