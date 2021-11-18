"""This contains functions to manipulate reaclib v2 data file"""

import struct
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

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

def load_pn(infile):
	"""
	Load pn table and write it to an array of dictionaries
	
	Parameters:
	   infile ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	datafrdmqrpa=[]
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1[0:len(line1)-1]
		val = line1.split()
		datafrdmqrpa.append({'Z': int(val[0]), 'N': int(val[1]), 'A': int(val[2]), 'P0n': float(val[3]), 'P1n': float(val[4]), 'P2n': float(val[5]), 'P3n': float(val[6]), 'P4n': float(val[7]), 'P5n': float(val[8]), 'P6n': float(val[9]), 'P7n': float(val[10]), 'P8n': float(val[11]), 'P9n': float(val[12]), 'P10n': float(val[13]), 'E_n': float(val[14]), 'n': float(val[15]), 'exp': int(val[16])})
	return datafrdmqrpa
def load_t12(infile):
	"""
	Load pn table and write it to an array of dictionaries
	
	Parameters:
	   infile ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	datafrdmqrpa=[]
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1[0:len(line1)-1]
		val = line1.split()
		datafrdmqrpa.append({'Z': int(val[0]), 'N': int(val[1]), 'A': int(val[0])+int(val[1]), 'T12': float(val[2])})
	return datafrdmqrpa

datafrdmqrpa_pxn = load_pn("pn-frdm2012-sdn-gtff-beoh350.dat")
datafrdmqrpa_t12 = load_t12("tlifminusff-beta-2018.dat")

datafrdmqrpa_pxn_t12 = []
for i in range(len(datafrdmqrpa_pxn)):
	match_entry = next((item for item in datafrdmqrpa_t12 if (item["Z"] == datafrdmqrpa_t12[i]["Z"] and item["N"] == datafrdmqrpa_t12[i]["N"])), None)
	if (match_entry==None):
		print("Error",datafrdmqrpa_pxn[i]["A"],datafrdmqrpa_pxn[i]["Z"]) 
	datafrdmqrpa_pxn_t12.append({'Z': datafrdmqrpa_pxn[i]['Z'], 'N': datafrdmqrpa_pxn[i]['N'], 'A': datafrdmqrpa_pxn[i]['A'], 'P0n': datafrdmqrpa_pxn[i]['P0n'], 'P1n': datafrdmqrpa_pxn[i]['P1n'], 'P2n': datafrdmqrpa_pxn[i]['P2n'], 'P3n': datafrdmqrpa_pxn[i]['P3n'], 'P4n': datafrdmqrpa_pxn[i]['P4n'], 'P5n': datafrdmqrpa_pxn[i]['P5n'], 'P6n': datafrdmqrpa_pxn[i]['P6n'], 'P7n': datafrdmqrpa_pxn[i]['P7n'], 'P8n': datafrdmqrpa_pxn[i]['P8n'], 'P9n': datafrdmqrpa_pxn[i]['P9n'], 'P10n': datafrdmqrpa_pxn[i]['P10n'], 'E_n': datafrdmqrpa_pxn[i]['E_n'], 'n': datafrdmqrpa_pxn[i]['n'], 'exp': datafrdmqrpa_pxn[i]['exp'],'T12':match_entry['T12']}) 

np.save("datafrdmqrpa_pxn_t12.npy",datafrdmqrpa_pxn_t12)

def drawbox(N,Z,fcolor='None',ecolor='gray', falpha = 1):
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

for i in range(len(datafrdmqrpa_pxn_t12)):
	plt.gca().add_patch(drawbox(datafrdmqrpa_pxn_t12[i]["N"],datafrdmqrpa_pxn_t12[i]["Z"],fcolor='gray',ecolor='None',falpha = 0.5))

plt.xlabel('Neutron number, $N$')
plt.ylabel('Proton number, $Z$')
plt.xlim([0,250])
plt.ylim([0,136])
plt.show()
