"""This contains functions to manipulate reaclib v2 data file"""

import struct
import numpy as np
import re

elements={"h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, 
"si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
 "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34, "br": 35, "kr": 36, "rb": 37, "sr": 38, "y": 39,
  "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44, "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
   "i": 53, "xe": 54, "cs": 55, "ba": 56, "la": 57, "ce": 58, "pr": 59, "nd": 60, "pm": 61, "sm": 62, "eu": 63, "gd": 64, "tb": 65,
    "dy": 66, "ho": 67, "er": 68, "tm": 69, "yb": 70, "lu": 71, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77, "pt": 78,
     "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85, "rn": 86, "fr": 87, "ra": 88, "ac": 89, "th": 90, "pa": 91,
      "u": 92, "np": 93, "pu": 94, "am": 95, "cm": 96, "bk": 97, "cf": 98, "es": 99, "fm": 100, "md": 101, "no": 102, "lr": 103, "rf": 104,
       "db": 105, "sg": 106, "bh": 107, "hs": 108, "mt": 109, "ds": 110, "rg": 111, "cn": 112, "nh": 113, "fl": 114, "mc": 115, "lv": 116, "ts": 117, "og": 118,
       "al-":119, "al*": 120}

def getnamebyz(z):
	"""
	Get element name by atomic number Z
	
	Parameters:
	   z ( int ): Atomic number Z
	"""
	keys=elements.keys()
	values=elements.values()
	return keys[values.index(z)]

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

def load_txt(infile):
	"""
	Load nubase file and write it to an array of dictionaries
	
	Parameters:
	   input ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	reaclib=[]
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1[0:len(line1)-1]
		if len(line1)<220:
			for i in range(220-len(line1)):
				line1 = line1+" "
		line1 = bytes(line1,encoding='utf8')

		(A,empty1,Zi,empty2,Ael,s_type,empty3,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,empty4,dT12,Jpi,Ensdf_year,Discov_year,BR) = struct.unpack("3s1s4s3s5s1s1s13s11s12s11s2s1s1s9s2s1s7s14s2s4s90s22x",line1)
		(A,empty1,Zi,empty2,Ael,s_type,empty3,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,empty4,dT12,Jpi,Ensdf_year,Discov_year,BR) = map(lambda x: x.decode('utf-8').strip(),(A,empty1,Zi,empty2,Ael,s_type,empty3,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,empty4,dT12,Jpi,Ensdf_year,Discov_year,BR))
		print(A,Ael,T12,dT12,BR)

load_txt('nubase_3.mas20.txt')