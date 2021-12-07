import numpy as np 
import pandas as pd

### Loop the data lines
with open("rae9-s1.data", 'r') as temp_f:
    # get No of columns in each line
    col_count = [ len(l.split()) for l in temp_f.readlines() ]
    #^no argment for split(), splits based on white space

### Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
column_names = range(0, max(col_count)) 
#^tells pandas what the max number of columns is 

### Read csv
df = pd.read_csv("rae9-s1.data", header=None, delimiter="\s+|;|:", names=column_names,engine="python")
cs=df.iloc[8]
print(cs)

class Input:
    
    # Constructor
    def __init__(self, filename):
        pass