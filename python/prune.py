# The purpose of this script is to take in a text file or multiple text files output by 
# NeutReflectometry, sift through its contents, and then output only neutron captures 
# and their child events.

# You should provide the filename as a command line argument, with at least two arguments. 
# The last will be the output.
# python3 prune.py input.txt output.txt
# python3 prune.py input1.txt intput2.txt input3.txt output.txt

# Import libraries -----------------------------------------------------
import sys
import os
import pandas as pd

# Take inputfile from cl and read it. ----------------------------------
# Start by initializing some stuff to edit.
infilenames = []
infilecontents = pd.DataFrame()
# Go through all the files.
for arg in sys.argv:
    if (arg != __file__ and arg != os.path.basename(__file__) and arg is not sys.argv[-1]): #Don't include the running of the script as an argument.
        infilenames.append(arg)                                 #note all our filenames
        infilecontents = pd.concat(                             #Combine all the data, but only the parts we need.
            [infilecontents,pd.read_csv(arg,delim_whitespace=True,usecols=['EV','TS','P','Type','E1','D3','time1','nCap'])]
            )

print(infilecontents)
#infilecontents = [x.replace('\n','').rstrip() for x in infilecontents] #Remove newlines at the end of each list entry.

#Add some header info before I modify the data. ------------------------
print("# Un-pruned file(s):", *infilenames, sep = ' ')
print("# Number of entries before pruning:",len(infilecontents.index))

#Here's the main the pruning process. ----------------------------------

outfilecontents = pd.DataFrame() #initialize this so I can add to it.
for value in infilecontents.query('nCap == 1')['EV']:
    outfilecontents = pd.concat([outfilecontents,                                                          #append something to the output contents.
        infilecontents[infilecontents.query('nCap == 1 and EV == @value').index[0]:].query('EV == @value')  #Include the whole event, starting from the initial capture
    ])

# pruning done, back to printing. --------------------------------------
print("# Capture events found:",len(infilecontents.query('nCap == 1').index)) #Print this before I change it
print("# Total lines output:",len(outfilecontents.index))
#print(open(sys.argv[1]).readlines()[1].replace('\n',''),'chainkey') #include the table header once.
outfilecontents.to_csv(sys.argv[-1])