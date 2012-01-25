####################################################################################################
###   Cluster.py :: Calculate Cluster and Aggregate data from LAMMPS data files output           ###
###   Written in Python ver 2.6 and above. *~* Clusterizer *~*                                   ###
###											         ###	
###								                                 ###
###   Funding gratefully provided by Unilever plc and University of Southampton.                 ###
###                                                                                              ###
###   All code is "as is" and there is no implied warranty.                                      ###
###                                                                                              ###
###   Brett Donovan (2012)	                                                                 ###
###                                                                                              ###
####################################################################################################

import cfg
import math
import time
import sys
import os
import re
from collections import defaultdict, deque
from array import array

def PullInSnapShot(solutemols, atomspersolmol):

	### We assume that solute molecules come first, we like to know how many atoms are in each molecule
	### so that we can finish when we have found all our molecules. Also id's should all be in order.
	
	atomtomoldict = {}
        atompositions = {}
	atomtomoltype = {}
	atompositions = {}
	timeclusterlist = []
	dim_vect = [] ## Store x/y of each frame as a list
	vol_vect = [] ## Store the volume of each frame in range [554, 545, 5345, 534, 34534....]
	linecnt = 0
	collecttime = collectdata = gotnumtime = collectbox = False
	f = open(cfg.inputtraj)
	for line in f:
		if "TIMESTEP" in line:
			collecttime = True
		if collecttime and not collectdata and not gotnumtime and re.search('\d+', line):
			time = int(re.search('\d+', line).group())
			gotnumtime = True
		if collecttime and "ITEM: ATOMS" in line:
			collectdata = True
		if collecttime and "ITEM: BOX BOUNDS" in line:
			collectbox = True
		if collectbox and re.search('\d+', line):
			if (linecnt == 0):
				x = abs(float(line.split()[0])*2.0)
			if (linecnt == 1):
				y = abs(float(line.split()[0])*2.0)
			if (linecnt == 2):
				z = abs(float(line.split()[0])*2.0)
			linecnt+=1
			#print "append", time, line, linecnt
		if (linecnt >= 3) and collectbox:
			#print "tot ", x, y, z
			vol = x*y*z
			linecnt = 0
			collectbox = False
		if collectdata and re.search('\d+', line):
			#print line
			#raw_input("")
			if (int(line.split()[2]) <= int(solutemols)):
				atom_data = line.split()
				atom = int(atom_data[0])
				atomtype = int(atom_data[1])
				atom_mod = modulus(atom, atomspersolmol)
				mol = int(atom_data[2])
				vectpos = atom_data[3], atom_data[4], atom_data[5]
				atomtomoldict[atom] = mol
        			atomtomoltype[atom] = atomtype
        			atompositions[atom] = vectpos		
				if int(line.split()[2]) == int(solutemols) and (atom_mod == atomspersolmol):
					collecttime = collectdata = gotnumtime = False

					if (time >= cfg.starttime) and (time <= cfg.finishtime):
 						dim_vect = [x, y, z]
						vol_vect.append(vol) 
						#print dim_vect, vol_vect
						print "Finished frame: ", time
						clusterdict = defaultdict(deque)
						LoopOverMols(clusterdict, atompositions, dim_vect)
						timeclusterlist.append(clusterdict)
						monomers, clustercount = CountMonomersClusters(clusterdict)
			
	#print len(vol_vect)			
	averagevolume = reduce(lambda a,b: a+b, vol_vect)/len(vol_vect)
	print "The average volume of the system: ", averagevolume
	Pn, elements = HistogramClusters(timeclusterlist, cfg.monomershistogram)
	#ComputeCMC(Pn, elements, averagevolume)
	print "vectors: ", dim_vect, "\n"
	f.close()		
	return atomtomoldict, atompositions, atomtomoltype, atompositions, dim_vect, time

def modulus(atom, atomspermol):
	op = int(atom-1) % (int(atomspermol))
	return op+1

def R2MinImageCo(atom1, atom2, dim_vect, atompositions):
	### Computes the R2 (radius squared) in the minimum image convention scheme.
	x = float(atompositions[atom1][0]) - float(atompositions[atom2][0])
	y = float(atompositions[atom1][1]) - float(atompositions[atom2][1])
 	z = float(atompositions[atom1][2]) - float(atompositions[atom2][2])
	x = x - dim_vect[0] * round(x/dim_vect[0])
	y = y - dim_vect[1] * round(y/dim_vect[1])
	z = z - dim_vect[2] * round(z/dim_vect[2])
	R2 = x**2 + y**2 + z**2
	return R2

def RR(mastermol, currentmol, atompositions, dim_vect, clusterdict):  
	atomcurrent = (currentmol-1)*int(cfg.atomspermol) + int(cfg.interestatom)
	for moltarget in range(1, int(cfg.nummols)+1):
		atomtarget = (moltarget-1)*int(cfg.atomspermol) + int(cfg.interestatom)  ### Compute the atom of interest on each molecule
		R2 = R2MinImageCo(atomcurrent, atomtarget, dim_vect, atompositions)
		R = math.sqrt(R2)
		if R < float(cfg.maxcutoff):
			#print atomcurrent, atomtarget, R
			#print dim_vect
			#raw_input("")
			success = AddToMaster(clusterdict, mastermol, moltarget)
			if success:
				RR(mastermol, moltarget, atompositions, dim_vect, clusterdict)

def isadup(clusterdict, value):
	
	### Search over all keys and values associated with every key to see whether we can

	duplicate = False
	if any(value in deq for deq in clusterdict.itervalues()):
		duplicate = True
	if value in clusterdict:
		duplicate = True
	return duplicate

def AddToMaster(clusterdict, master_mol, mol_id):

	### We use the initial molecule as a designator and build a list of molecules which are 
        ### on the same cluster. We can then scan through the list of these.
	### Format: Monomers: 1:1 2:2 99:99, Clusters: 1:[1, 2, 4, 5, 9, 11]
	### For ease of reading and counting we do duplicate the key in the value dictionary pair (means
	### we get the monomer count easily).

	if not isadup(clusterdict, mol_id):
		clusterdict[master_mol].append(mol_id)
		return True
	else:
		return False

def CountMaxClusterSize(clusterdict):
	maxsize = 0
	for deq in clusterdict.itervalues():
		if (len(deq) > maxsize):
			maxsize = len(deq)
	return maxsize

def CountMonomersClusters(clusterdict):
	monomercount = clustercount = 0
	for deq in clusterdict.itervalues():
		if len(deq) == 1:
			monomercount+=1
			#print "Monomers found! Yippee!"
		else:
			clustercount+=1
			#print "Clusters found! Oh Yes!"
	return monomercount, clustercount



def HistogramClusters(clusterlist, monomers):

	### Monomers is a boolean and switches on the results for monomers
	
	summation = 0.0
	endbuffer = 0
	f = open(cfg.outputhist, "w")
	maxclustersize = 0
	for clusterdict in clusterlist:
		maxclustersizedict = CountMaxClusterSize(clusterdict)
		if (maxclustersizedict > maxclustersize):
			maxclustersize = maxclustersizedict
	histarray = []
	for agg in range(maxclustersize+1):
		histarray.append(0)
	for clusterdict in clusterlist:
		for deq in clusterdict.itervalues():
			length = len(deq)
			#print length
			histarray[length]+=1
	for agg in range(maxclustersize+endbuffer+1):
		summation += histarray[agg]
	## Calculate P(n)
	for agg in range(maxclustersize+endbuffer+1):
		histarray[agg] = histarray[agg]/summation

	for agg in range(1, maxclustersize+endbuffer+1):
		if not monomers:
			if (agg > 1):
				line = "%d %f\n" % (agg, histarray[agg])
				f.write(line)
		else:
			line = "%d %f\n" % (agg, histarray[agg])
			f.write(line)
	f.close()
	return histarray, (maxclustersize+1) ## Return P(n)
	

def ComputeCMC(Pn_array, elements, averagevolume):

	### Figure out what the CMC based on integration of P(n) from n=1 to Peak
	### Compute Peaks - assume there is more than one and pick the highest one.
	
	peakdict = {}
	for agg in range(2, elements):
		Pn = Pn_array[agg]
		if (Pn_array[agg-1] < Pn) and (Pn_array[agg+1] < Pn):
			peakdict[agg] = Pn_array[agg]

	lastpeakagg = peakdict.items()[-1][0]

	### Lets integrate this up
	### bin width = 1

	aggtotal = 0
	for intrange in range(2, lastpeakagg):
		aggtotal += Pn_array[intrange]
		

	### Lets compute the volume fraction
	### Approximate volume of molecules (cfg.thrvolumemol)	
	### Compute average of total volume over all frames

	cmcvolfrac = (aggtotal*float(cfg.thrvolumemol))/float(averagevolume)
	averagevolinm3 = averagevolume * 1E-30
	molconcentration = aggtotal / (6.022E23)

	print "Integral under Pn(n) is: ", aggtotal
	print "Volume fraction at cmc v/v: ", cmcvolfrac 
	print "CMC in trad. units: ", molconcentration/averagevolinm3, "mM"
		
	
def LoopOverMols(clusterdict, atompositions, dim_vect):
	for mol in range(1, int(cfg.nummols)+1):
		if not isadup(clusterdict, mol):
			RR(mol, mol, atompositions, dim_vect, clusterdict)

	### Do our analysis on clusters
	### Compute number of monomers

		
def main():
	print "Clustering algorithm...."
	#clusterdict = defaultdict(deque)
	atomtomoldict, atompositions, atomtomoltype, atompositions, dim_vect, time = PullInSnapShot(cfg.nummols, cfg.atomspermol)
	#LoopOverMols(clusterdict, atompositions, dim_vect)
	#CountMonomers(clusterdict)
	#print clusterdict


if __name__ == "__main__":
	main()


