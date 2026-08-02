# clean.py
''' 
This script will import and clean the positional data from four ground stations along the 
cascades region.
Jun 20, 2026
Version 2
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



