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

def plot_complement_data():
	magic_num = [2, 8, 20, 28, 50, 82, 126]
	for i in magic_num:
		plt.axhline(y=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axhline(y=i-0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i-0.5,color='b',linestyle='--',linewidth=0.2)
	data_bound = np.load("data_bound.npy",allow_pickle='TRUE')
	for i in range(len(data_bound)):
		plt.gca().add_patch(drawbox(data_bound[i]["N"],data_bound[i]["Z"],fcolor='gray',ecolor='None',falpha = 0.5))
	nubase_stable_add_to_iaea_crp = np.load("nubase_stable_add_to_iaea_crp.npy",allow_pickle='TRUE')
	for i in range(len(nubase_stable_add_to_iaea_crp)):
		plt.gca().add_patch(drawbox(nubase_stable_add_to_iaea_crp[i]["N"],nubase_stable_add_to_iaea_crp[i]["Z"],fcolor='k',ecolor='None',falpha = 1))
	datafrdmqrpa_pxn_t12_add_to_iaeacrp_bdn_sep = np.load("datafrdmqrpa_pxn_t12_add_to_iaeacrp_bdn_sep.npy",allow_pickle='TRUE')
	for i in range(len(datafrdmqrpa_pxn_t12_add_to_iaeacrp_bdn_sep)):
		plt.gca().add_patch(drawbox(datafrdmqrpa_pxn_t12_add_to_iaeacrp_bdn_sep[i]["N"],datafrdmqrpa_pxn_t12_add_to_iaeacrp_bdn_sep[i]["Z"],fcolor='g',ecolor='k',falpha = 1,linewidth=0.001))
	nubase_bminus_add_to_iaeacrp_bdn = np.load("nubase_bminus_add_to_iaeacrp_bdn.npy",allow_pickle='TRUE')
	for i in range(len(nubase_bminus_add_to_iaeacrp_bdn)):
		plt.gca().add_patch(drawbox(nubase_bminus_add_to_iaeacrp_bdn[i]["N"],nubase_bminus_add_to_iaeacrp_bdn[i]["Z"],fcolor='r',ecolor='k',falpha = 1,linewidth=0.001))
	iaea_crp_bdn = np.load("iaea_crp_bdn.npy",allow_pickle='TRUE')
	for i in range(len(iaea_crp_bdn)):
		plt.gca().add_patch(drawbox(iaea_crp_bdn[i]["N"],iaea_crp_bdn[i]["Z"],fcolor='y',ecolor='k',falpha = 1,linewidth=0.001))
	plt.xlabel('Neutron number, $N$')
	plt.ylabel('Proton number, $Z$')
	plt.xlim([9.5,200])
	plt.ylim([9.5,116])

plot_complement_data()
plt.show()