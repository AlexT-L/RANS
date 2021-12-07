import numpy as np 
import pandas as pd

###Get max number of columns in a row
def max_no_cols(file):
    ### Loop the data lines
    with open(file, 'r') as temp_f:
        # get No of columns in each line
        col_count = [ len(l.split()) for l in temp_f.readlines() ]
        #^no argment for split(), splits based on white space
    return max(col_count)


### Read .data file

def read(file,max_cols):
    df = pd.read_csv(file, header=None, delimiter="\s+|;|:", names=range(max_cols),engine="python")
    return df

#class Input:
    
    # Constructor
    #def __init__(self, filename):
     #   pass
max=max_no_cols("rae9-s1.data")
print (read("rae9-s1.data",max))