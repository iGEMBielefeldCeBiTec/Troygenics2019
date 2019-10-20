### iGEM Bielefeld-CeBiTec 2019, Isabel Conze ###
### iconze@cebitec.uni-bielefeld.de ###
### v1 ###

__usage__ = """
		Calculate the time it takes until Saccharomyces serevisae gains resistance against a Cas13a system
		"""
#imports:
import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt

#mutation frequencies: S. cerevisiae
sbm=1.67*10**(-12)
patta=sbm*0.063
pcggc=sbm*0.152
pcagt=sbm*0.182
pactg=sbm*0.11
pagtc=sbm*0.144
pctga=sbm*0.35
indel_all=5.03*10**-12 #probability: insertion or deletion
genic=0.75 #percentage genic area
gen_indel=0.53 #percentage of indels found in genic areas
indel=int(indel_all*genic*gen_indel) #likelihood of indel taking place in a genic area


mwa=patta+pactg+pagtc #mutation frequency "away from A"
mwt=patta+pactg+pagtc #mutation frequency "away from T"
mwc=pcggc+pcagt+pctga #mutation frequency "away from C"
mwg=pcggc+pcagt+pctga #mutation frequency "away from G"
mwga=mwg+mwc
mwag=mwt+mwa

#relevant for mutation probability:
gs=36456735.001  #genome size (bp)
gl=21 #gRNA length
mn=1.0 #mutations needed for resistance against Cas13a

gc01=0.476 #gc content gRNA1
gc02=0.476 #gc content gRNA2
gc03=0.571 #gc content gRNA3
gc04=0.38 #gc content gRNA4
gc05=0.476 #gc content gRNA5
gc06=0.524 #gc content gRNA6
gc07=0.476 #gc content gRNA7
gc08=0.5 #gc content gRNA8
gc09=0.5 #gc content gRNA9
gc10=0.5 #gc content gRNA10

p_lalo=6.3*10**(-8) #likelihood of large losses (>1 kb)
lenBI=5844 #length Cas and gRNAs with their respective promoters and terminators
lalo=p_lalo*(lenBI/gs)

#mutation frequency for each gRNA: (GC-relations for 7 possible chosen gRNAs)
mr01=((gc01*mwga)+((1-gc01)*mwag)+indel)*(gl/gs)+lalo
mr02=((gc02*mwga)+((1-gc02)*mwag)+indel)*(gl/gs)+lalo
mr03=((gc03*mwga)+((1-gc03)*mwag)+indel)*(gl/gs)+lalo
mr04=((gc04*mwga)+((1-gc04)*mwag)+indel)*(gl/gs)+lalo
mr05=((gc05*mwga)+((1-gc05)*mwag)+indel)*(gl/gs)+lalo
mr06=((gc06*mwga)+((1-gc06)*mwag)+indel)*(gl/gs)+lalo
mr07=((gc07*mwga)+((1-gc07)*mwag)+indel)*(gl/gs)+lalo
mr08=((gc08*mwga)+((1-gc08)*mwag)+indel)*(gl/gs)+lalo
mr09=((gc09*mwga)+((1-gc09)*mwag)+indel)*(gl/gs)+lalo
mr10=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr11=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr12=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr13=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr14=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr15=((gc10*mwga)+((1-gc10)*mwag)+indel)*(gl/gs)+lalo
mr=[mr01,mr02,mr03,mr04,mr05,mr06,mr07,mr08,mr09,mr10,mr11,mr12,mr13,mr14,mr15]
print mr

#situation in the beginning: 
g=1 #generation
startOD=0.4
cs=0.4*10**7 #no of cells: sensitive
cr=0 #no of cells: resistant
c1=0.4*10**7 #cell count, beginning
K=15*10**7 #maximal cell number
K1=15*10**7
gt=2.3544 #generation time (hours)
t=0 #time (hours since first infection)
ts=0 #time (hours since first resistant cell)
r=1.06 #growth rate without Cas13a
r1=0.12 #growth rate with Cas13a
A=(K-c1)/c1 #for growth curve

#starting -lists:
days=[]
cells=[]
cellss=[]
cellsr=[]
days.append(t)
cellss.append(cs)
cellsr.append(cr)
pop=cs+cr
cells.append(pop)
i=1
a=0
crr=0
p_list=[]
mri=[]
trt=[]

#look at the mutation frequency of each gRNA:
for i in range(0,15): #i=number of gRNAs
	mri.append(mr[i]) #list of all "active" gRNAs
	print mri
	p=np.prod(mri)
	a=1.0
	b=2
	cs=0
	cells=[]
	cellss=[]
	cellsr=[]
	days=[]
	t=0
	cr=0
	K=15*10**7
	p_list.append(p)
	print p, "p", gt, "gt", cs, "cs"
	
	while cr==0: #while there is no resistant (=1) cell in the population, continue
		t+=gt
		cs=int(c1*(np.exp(np.log(K/c1)*(1-np.exp(-r1*t)))))
		g+=1
		days.append(t)
		cellss.append(cs)
		
		if cs*p >= 0.85: #the likelihood of gaining resistance is larger than 85%
			gr=g
			tr=t
			cr=1
			c2=1
			cellsr.append(cr)
			cells.append(cr+cs)
			trt.append(t)
			ts+=gt
			t+=gt
			
			break
			
		elif t > 26280:
			print "It takes more than 3 years to gain resistance!"
			cells.append(cs+cr)
			cellsr.append(cr)
			gr="more than 2000 h"
			tr="more than" + str(g)
			cr=1
		
		elif cs > (K-1000): #No resistant cells up to stationairy phase, wash and resuspend in same amount of fresh medium!"
			b+=1
			K=2*K
			cellsr.append(cr)
			cells.append(cs+cr)
			
		else: 
			cellsr.append(cr)
			cells.append(cs+cr)
		
	fig=plt.figure()
	plt.plot(days, cellss)
	plt.title("Growth Curve S. cerevisiae (Modeled)")
	plt.xlabel("time [hours]")
	plt.ylabel("number of cells")
	plt.show()
	plt.savefig(str(i)+"growth_SC")
	plt.close()

		
	g=1 #generation
	cs=10 #no of cells: sensitive
	cr=0 #no of cells: resistant
	pop=cr+cs
	c1=10 #cell count, beginning
	t=0
	ts=0
	days=[]
	cells=[]
	cellss=[]
	cellsr=[]
	days.append(t)
	cellss.append(cs)
	cellsr.append(cr)
	cells.append(pop)
	pops=[0] #population: sensitive cells
	popr=[0] #population: resistant cells
	rg=[]
	a+=1
	b=0
	crr=0
	cs1=0
	cr1=0
	gr="nd"
	tr="nd"
	

gS=range(1,len(p_list)+1)

years = gS
visitors = p_list
index = np.arange(len(gS))
bar_width = 0.9
plt.bar(index, p_list, bar_width,  color="deepskyblue")
plt.yscale("log")
plt.xlabel("Number of gRNAs")
plt.ylabel("Probability of gaining resistance")
plt.show()
plt.savefig("Probability_Gom.png")
plt.close()

plt.plot(gS, trt, color="deepskyblue")
plt.xlabel("Number of gRNAs")
plt.ylabel("Time until first resistant cell arises [h]")
plt.savefig("TimeSC_Gom.png")
plt.show()
