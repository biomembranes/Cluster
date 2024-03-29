import math

### Clustering inputs

inputtraj = "dump.trj"
nummols = 1280
atomspermol = 9
interestatom = 5
maxcutoff = 12 ## in A
monomershistogram = True  ## Display monomer count in histograms
thrvolumemol = 600 ## theoretical volume of the molecule in A^3

### Time range

starttime  = 0
finishtime = 2000000


### Histogram Output

outputhist = "hist.dat"

### Data Formats

ROWDATA  = 7 ## This is the number of lumps of data we take from LAMMPS data file
ROWBOND  = 4 ## Bond data from the LAMMPS data file
ROWANGLE = 5 ## Angle data from the LAMMPS data file
