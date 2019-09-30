import pandas as pd
import numpy as np
import sqlite3

def create_sql_database():

    pd_data = pd.read_csv("KeckSpreadsheet.csv")
    index = np.arange(len(pd_data.index))
    pd_data["index"] = index

    conn = sqlite3.connect("keckdatabase.db")
    c = conn.cursor()
    
    # Object
    # id, name, ra, dec, champion
    object_headers = ["index", "HIP", "RA (deg)", "dec (deg)", "Champion"]
    object_values = pd_data[object_headers].values.tolist()
    c.executemany('INSERT INTO Object VALUES (?, ?, ?, ?, ?)', object_values)

    # Measurements
    # id, type, distance, gaiadistance, sptype, propermotion, vmag, bminusv
    measurement_headers = ["index", "Target Type", "dist", "Gaia DR2 dist", "SpType", "PM (mas/yr)", "V", "(B-V)"]
    measurement_values = pd_data[measurement_headers].values.tolist()
    c.executemany('INSERT INTO Measurements VALUES (?, ?, ?, ?, ?, ?, ?, ?)', measurement_values)

    # Observations
    # id, followup, candidates, reductionnotes, comment
    observation_headers = ["index", "Need followup?", "Candidates?", "Reduction Notes", "Notes"]
    observation_values = pd_data[observation_headers].values.tolist()
    c.executemany('INSERT INTO Observations VALUES (?, ?, ?, ?, ?)', observation_values)

    conn.commit()
    conn.close()
    
        
if __name__ == "__main__":
    create_sql_database()