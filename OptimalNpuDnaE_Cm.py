### iGEM Bielefeld-CeBiTec 2019, Isabel Conze ###
### iconze@cebitec.uni-bielefeld.de ###
### v1 ###

__usage__ = """
		Find the ideal split point for the splicing of Chloramphenicolacetyltransferase with Npu DnaE
		"""

#imports:
import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt

#Amino Acid Sequence
Seq=("MEKKITGYTTVDISQWHRKEHFEAFQSVAQCTYNQTVQLDITAFLKTVKKNKHKFYPAFIHILARLMNAHPEFRMAMKDGELVIWDSVHPCYTVFHEQTETFSSLWSEYHDDFRQFLHIYSQDVACYGENLAYFPKGFIENMFFVSANPWVSFTSFDLNVANMDNFFAPVFTMGKYYTQGDKVLMPLAIQVHHAVCDGFHVGRMLNELQQYCDEWQGGA")
Len=range(1,len(Seq))

#Sequence of Structures
struc=["DISQW","RKEH","HFEAFQ","ITAFLKTVKKNK","KFYPAFIHILARLMN","HPEF","DFRQFLHIYSQDVACYG","DGFHVGRMLNELQQYCDEW","YTT","NQTVQLD","MAMK","ELVIW","HPCYTVFH","TFSSLW","MFFVSA","SFDLNV","VFTM","YTQ", "KVLMPLAIQVH"]

#sequence of structure to index
structure=[]
for i in range (0,len(struc)):
	if struc[i] in Seq:
		structure.append(range(Seq.index(struc[i]),Seq.index(struc[i])+len(struc[i])))
print structure, "structure!"

proa=[]
pro_com=[]
pro_comb=[]

#All possible "cut" sequences: pro_com
for i in range(0,len(Seq)-6):
	proa=[Seq[i],Seq[i+1], Seq[i+2], Seq[i+3], Seq[i+4], Seq[i+5]]
	print proa
	pro_com.append(proa)
	pro_comb.append(proa)
	


#how good is the amino acid in that position
D=[1.39, 1.52, 0, 0, 0, 1.26]
E=[1.26,3.51,0.07, 0, 0, 5.5]
N=[0.93,1.26,2.19, 0,0,8.61]
Q=[1.39,0.86,0.33, 0,0,0.27]
H=[1.06,1.26,2.25, 0,0,1.59]
K=[1.59,0.93,4.44, 0,0.07,0]
R=[1.81,0.2,0.86, 0,0,0]
S=[0.99,1.24,1.37, 0,0,0.07]
C=[0.46,0.4,0.8,50,0.07,0.07]
T=[0.63,0.89,1.59, 0,0,0.96]
P=[0.8,1.52,0,0,0,0]
G=[4.67,2.68,0.17,0,0,0.03]
A=[0.63,0.83,1.92,0,0,0.23]
V=[0.56,0.5,0.23,0,0.03,0.1]
I=[0.07,0.53,0,0,0,0.4]
L=[0.11,0.57,0.75,0,0,0.82]
M=[0.13,1.13,0.86,20,2.19,1.46]
F=[0.13,0.73,1.26,0,0.07,2.25]
Y=[0.2,0.6,2.05,0,0.13,5.23]
W=[0.07,0.4,0.99,0,29.42,0.07]

#all amino acids
AAs=["D","E","N","Q","H","K","R","S","C","T","P","G","A","V","I","L","M","F","Y","W"]
AA2=[D,E,N,Q,H,K,R,S,C,T,P,G,A,V,I,L,M,F,Y,W]

#empty list for filling with "how good is this combination"
pro_num=[[0,0,0,0,0,0] for x in xrange(0,len(pro_com))]

#add: "how good is this amino acid in that combination" (final list: pro_num)
for i in range (0,len(pro_com)):
	for x in range(0,len(AAs)):
		if AAs[x] in pro_com[i]:
			c=AAs[x]
			a=[e for e, g in enumerate(pro_com[i]) if g == c] 
			for z in range(0,len(a)):
				f=a[z]
				rn=AA2[x][f] 
				pro_num[i][a[z]]=rn
				

splitval=[]
splitvalb=[]

#how good is each combination (final list: splitval/splitvalb)
for i in range(0,len(pro_num)):
	d=sum(pro_num[i])
	splitval.append(d)
	splitvalb.append(d)


#sort the combination-sums (final list: sumsort)
sumsort=sorted(splitval,reverse=True)


indexes=[]
indb=[]
#gives a list of indexes in order (e.g. the second combination is the best, adds 1 as the first element of the list)
for i in range(0,len(sumsort)):
	h=sumsort[i]
	imax=splitvalb.index(h) 
	indexes.append(imax)
	indb.append(imax)
	splitvalb[imax]=0

print "indexes", indexes
sort=[]
#sorted list of amino-acid-combinations
for i in range (0, len(indexes)):
	j=indexes[i]
	sort.append(pro_com[j])
print sort

a=0	
#checks if the first amino acid is part of a relevant structure
while a != len(indexes):
	for i in range (0, len(structure)):
		if indexes[a] in structure[i]:
			indb[a]=0

	a+=1
l=len(indb)	
i=0
while(i<l):
	if(indb[i]==0):
		indb.remove (indb[i])	
		l-=1   
		continue
	i+=1
		
print "The best positions for intein splicing are", indb

fAAs=[]
for i in range (0,len(indb)):
	g=indb[i]
	a=pro_com[g]
	fAAs.append(a)
	
print "The AA-combinations of those are", fAAs

fobj = open("AminoAcids_Cm_withoutInsertion.txt", "w")
with open("AminoAcids_Cm_withoutInsertion.txt", "w") as fh:
	fh.write(str(fAAs))
	fh.write(str(Seq))
