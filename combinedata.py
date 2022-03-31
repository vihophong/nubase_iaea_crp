import struct
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

time_factor = {'s':1., 'y':31536000., 'ms': 0.001, 'd' : 86400., 'ky' : 31536000000, 'm' : 60., 'h': 3600.}

elements={"h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, 
"si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
 "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34, "br": 35, "kr": 36, "rb": 37, "sr": 38, "y": 39,
  "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44, "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
   "i": 53, "xe": 54, "cs": 55, "ba": 56, "la": 57, "ce": 58, "pr": 59, "nd": 60, "pm": 61, "sm": 62, "eu": 63, "gd": 64, "tb": 65,
    "dy": 66, "ho": 67, "er": 68, "tm": 69, "yb": 70, "lu": 71, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77, "pt": 78,
     "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85, "rn": 86, "fr": 87, "ra": 88, "ac": 89, "th": 90, "pa": 91,
      "u": 92, "np": 93, "pu": 94, "am": 95, "cm": 96, "bk": 97, "cf": 98, "es": 99, "fm": 100, "md": 101, "no": 102, "lr": 103, "rf": 104,
       "db": 105, "sg": 106, "bh": 107, "hs": 108, "mt": 109, "ds": 110, "rg": 111, "cn": 112, "nh": 113, "fl": 114, "mc": 115, "lv": 116, "ts": 117, "og": 118,
	   "119": 119,"120": 120,"121": 121,"122": 122,"123": 123,"124": 124,"125": 125,"126": 126,"127": 127,"128": 128,"129": 129,"130": 130,
	   "131": 131,"132": 132,"133": 133,"134": 134,"135": 135,"136": 136}

Zele = []
Zele.append("n")
for key in elements:
	Zele.append(key)

def getnamebyz(z):
	"""
	Get element name by atomic number Z
	
	Parameters:
	   z ( int ): Atomic number Z
	"""
	return Zele[z]

def getZ(input):
	"""
	Get atomic number Z by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -8888
	else:
		sep=re.split('(\d+)',input)
		if len(sep)==1:
			if sep[0]=="n":
				return int(0)
			elif (sep[0]=="p" or sep[0]=="d" or sep[0]=="t"):
				return int(1)			
			else:
				print("Something wrong! ",input)
		else:
			return int(elements[sep[0]])

def getA(input):
	"""
	Get mass number A by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -9999
	else:
		sep=re.split('(\d+)',input)
		if len(sep)==1:
			if sep[0]=="n":
				return 1
			elif sep[0]=="p":
				return 1
			elif sep[0]=="d":
				return 2
			elif sep[0]=="t":
				return 3
			else:
				print("Something wrong! ",input)
		else:
			return int(sep[1])

def drawbox(N,Z,fcolor='None',ecolor='gray', falpha = 1,linewidth=1):
	"""
	Draw box
	
	Parameters:
	   N ( int ): Neutron number
	   P ( int ): Proton number
	   ecolor ( str ): Color code
	"""
	rec = plt.Rectangle((N-0.5,Z-0.5),1,1,facecolor=fcolor,edgecolor=ecolor,alpha = falpha)
	#plt.text(N-0.4,Z-0.1,'$\mathregular{^{'+str(Z+N)+'}'+elements[Z]+'}$')
	return rec 

def load_iaea_crp(infile):
	"""
	Load IAEA_CRP file and write it to an array of dictionaries
	
	Parameters:
	   infile ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	iaea_crp_bdn = []
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1.strip()
		val = line1.split()
		# for i in val:
		# 	print(val)
		Z = float(val[1])
		A = float(val[2])
		N = A-Z
		liso = float(val[3])
		# Skipping isomer
		if (liso!=0):
			continue

		T12 = float(val[16])
		dT12 = float(val[17])
		P1n = float(val[18])
		dP1n = float(val[19])
		P2n = float(val[20])
		dP2n = float(val[21])
		P3n = float(val[22])
		dP3n = float(val[23])
		dT12hi = float(val[30])
		dP1nhi = float(val[31])
		dP2nhi = float(val[32])
		iaea_crp_bdn.append({"A":A,"Z":Z,"N":N, "T12":T12, "dT12":dT12, "dT12hi":dT12hi, "P1n": P1n, "dP1n": dP1n, "dP1nhi":dP1nhi, "P2n": P2n, "dP2n": dP2n, "dP2nhi":dP2nhi, "P3n": P3n, "dP3n": dP3n, "source":"iaeacrp"})
	np.save("iaea_crp_bdn_220327.npy",iaea_crp_bdn)
	print(iaea_crp_bdn)

load_iaea_crp('220327_listofeval_exp.txt')

def combinedata():
	iaea_crp_bdn = np.load("iaea_crp_bdn_220327.npy",allow_pickle='TRUE')
	# data_bound = np.load("data_bound.npy",allow_pickle='TRUE')
	nubase_stable = np.load("nubase_stable.npy",allow_pickle='TRUE')
	nubase_bminus = np.load("nubase_bminus.npy",allow_pickle='TRUE')
	#data not overlap with iaea_crp_bdn
	nubase_bminus_add_to_iaeacrp_bdn = []
	for i in range(len(nubase_bminus)):
		match_entry_nubase_bminus = next((item for item in iaea_crp_bdn if (item["Z"] == nubase_bminus[i]["Z"] and item["N"] == nubase_bminus[i]["N"])), None)
		if (match_entry_nubase_bminus==None):
			nubase_bminus_add_to_iaeacrp_bdn.append(nubase_bminus[i])
	np.save("nubase_bminus_add_to_iaeacrp_bdn_220327.npy",nubase_bminus_add_to_iaeacrp_bdn)
	nubase_stable_add_to_iaeacrp_bdn = []
	for i in range(len(nubase_stable)):
		match_entry_nubase_stable = next((item for item in iaea_crp_bdn if (item["Z"] == nubase_stable[i]["Z"] and item["N"] == nubase_stable[i]["N"])), None)
		if (match_entry_nubase_stable==None):
			nubase_stable_add_to_iaeacrp_bdn.append(nubase_stable[i])
	np.save("nubase_stable_add_to_iaeacrp_bdn_220327.npy",nubase_stable_add_to_iaeacrp_bdn)
	iaea_crp_nubase_combined = []
	for i in nubase_bminus_add_to_iaeacrp_bdn:
		iaea_crp_nubase_combined.append({"A":i["A"],"Z":i["Z"],"N":i["N"], "T12":i["T12"], "dT12":i["dT12"], "dT12hi":i["dT12"], "P1n": i["P1n"], "dP1n": i["dP1n"], "dP1nhi":i["dP1n"], "P2n": i["P2n"], "dP2n": i["dP2n"], "dP2nhi":i["dP2n"], "P3n": 0., "dP3n": 0., "source":"nubase"})
	for i in iaea_crp_bdn:
		iaea_crp_nubase_combined.append(i)
	np.save("iaea_crp_nubase_combined_220327.npy",iaea_crp_nubase_combined)
combinedata()

def plotcombineddata():
	magic_num = [2, 8, 20, 28, 50, 82, 126]
	for i in magic_num:
		plt.axhline(y=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axhline(y=i-0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i-0.5,color='b',linestyle='--',linewidth=0.2)
	nubase_stable_add_to_iaea_crp = np.load("nubase_stable_add_to_iaeacrp_bdn_220327.npy",allow_pickle='TRUE')
	# print(nubase_stable_add_to_iaea_crp)
	for i in range(len(nubase_stable_add_to_iaea_crp)):
		plt.gca().add_patch(drawbox(nubase_stable_add_to_iaea_crp[i]["N"],nubase_stable_add_to_iaea_crp[i]["Z"],fcolor='k',ecolor='None',falpha = 1))
	iaea_crp_nubase_combined = np.load("iaea_crp_nubase_combined_220327.npy",allow_pickle='TRUE')
	for i in range(len(iaea_crp_nubase_combined)):
		if (iaea_crp_nubase_combined[i]["source"]=="nubase"):
			plt.gca().add_patch(drawbox(iaea_crp_nubase_combined[i]["N"],iaea_crp_nubase_combined[i]["Z"],fcolor='r',ecolor='k',falpha = 1,linewidth=0.001))
		else:
			plt.gca().add_patch(drawbox(iaea_crp_nubase_combined[i]["N"],iaea_crp_nubase_combined[i]["Z"],fcolor='y',ecolor='k',falpha = 1,linewidth=0.001))
			plt.text(float(iaea_crp_nubase_combined[i]["N"])-0.2,float(iaea_crp_nubase_combined[i]["Z"])-0.2, str(int(iaea_crp_nubase_combined[i]["A"]))+getnamebyz(int(iaea_crp_nubase_combined[i]["Z"])).capitalize(),fontsize='xx-small')
	plt.xlabel('Neutron number, $N$')
	plt.ylabel('Proton number, $Z$')
	plt.xlim([9.5,200])
	plt.ylim([9.5,116])

plotcombineddata()
plt.show()
