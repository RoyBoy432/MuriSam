# -*- coding: utf-8 -*-
"""
Created on 2020/02/19

@author: Roy Moger-Reischer
"""

#%%
from __future__ import division
import re, os,sys, math, operator,random, copy,collections,time; import numpy as np; import pandas as pd; import csv
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
#%%
import vcf
#%%
#mytest=pd.read_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF')

mytest = vcf.Reader(open(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF','r'))
for r in mytest:
    print(r)
    

#%%
def count_fixed_mutations_3B(infile="",samplename=""):
    sample = vcf.Reader(open(infile,'r'))
    basis=dict()
    basis['in']=0;basis['del']=0;basis["SV"]=0;basis["AT_CG"]=0;basis["AT_GC"]=0;basis["AT_TA"]=0;basis["CG_GC"]=0;basis["CG_TA"]=0;basis["GC_TA"]=0
    for m in sample:
        if m.POS >= 76900 and m.POS <= 76909:#there is an oligo-T deletion in the MA_3B ancestor, but breseq doesn't always score it exactly the same, and therefore it wasn't subtracted out of the lines---hence I am subtracting it here
            continue
        else:
            if m.INFO["AF"][0]!=1.0:
                #print("not fixed")
                continue
            else:
                #print(m.ALT[0])
                if len(m.ALT[0])!=1 or len(m.REF)!=1:
                    #figure out if it's an in, a del, or a SV
                    if len(m.ALT[0])>100 or len(m.REF)>100:
                        basis["SV"]+=1
                    elif len(m.ALT[0]) > len(m.REF):
                        basis['in']+=1
                    elif len(m.ALT[0]) < len(m.REF):
                        basis['del']+=1
                    else:
                        continue
                else:
                    if m.REF == "A":
                        if m.ALT[0] == "C":
                            basis["AT_CG"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_GC"]+=1
                        else:
                            basis["AT_TA"]+=1
                    elif m.REF == "C":
                        if m.ALT[0] == "A":
                            basis["GC_TA"]+=1
                        elif m.ALT[0] == "G":
                            basis["CG_GC"]+=1
                        else:
                            basis["CG_TA"]+=1
                    elif m.REF == "G":
                        if m.ALT[0] == "A":
                            basis["CG_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["CG_GC"]+=1
                        elif m.ALT[0] == "T":
                            basis["GC_TA"]+=1
                    else:
                        if m.ALT[0] == "A":
                            basis["AT_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["AT_GC"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_CG"]
                    #figure out which of the 6 possible substs it was
            
    return basis

#%%
count_fixed_mutations_3B(r'C:\Users\rmoge\GitHub\MuriSam\MA_testing\3B_1_clean.VCF')

#%%
wholedict=dict()
for n in range(1,97):
    stringname="3B_" + str(n)
    thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_VCF\3B_' + str(n) + r'_clean.VCF'
    wholedict[stringname]=count_fixed_mutations_3B(thisfile)
MA_3B_df=pd.DataFrame.from_dict(wholedict,orient="index")
MA_3B_df.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_3B_VCF\MA_3B_summary.csv')


#%%
def count_fixed_mutations_s1(infile="",samplename=""):
    sample = vcf.Reader(open(infile,'r'))
    basis=dict()
    basis['in']=0;basis['del']=0;basis["SV"]=0;basis["AT_CG"]=0;basis["AT_GC"]=0;basis["AT_TA"]=0;basis["CG_GC"]=0;basis["CG_TA"]=0;basis["GC_TA"]=0
    for m in sample:
        if m.POS == 917603:#there is mutation at this position in the ancestor that was not caught by gdtools subtract because breseq didnt catch it
            continue
        else:
            if m.INFO["AF"][0]!=1.0:
                #print("not fixed")
                continue
            else:
                #print(m.ALT[0])
                if len(m.ALT[0])!=1 or len(m.REF)!=1:
                    #figure out if it's an in, a del, or a SV
                    if len(m.ALT[0])>100 or len(m.REF)>100:
                        basis["SV"]+=1
                    elif len(m.ALT[0]) > len(m.REF):
                        basis['in']+=1
                    elif len(m.ALT[0]) < len(m.REF):
                        basis['del']+=1
                    else:
                        continue
                else:
                    if m.REF == "A":
                        if m.ALT[0] == "C":
                            basis["AT_CG"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_GC"]+=1
                        else:
                            basis["AT_TA"]+=1
                    elif m.REF == "C":
                        if m.ALT[0] == "A":
                            basis["GC_TA"]+=1
                        elif m.ALT[0] == "G":
                            basis["CG_GC"]+=1
                        else:
                            basis["CG_TA"]+=1
                    elif m.REF == "G":
                        if m.ALT[0] == "A":
                            basis["CG_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["CG_GC"]+=1
                        elif m.ALT[0] == "T":
                            basis["GC_TA"]+=1
                    else:
                        if m.ALT[0] == "A":
                            basis["AT_TA"]+=1
                        elif m.ALT[0] == "C":
                            basis["AT_GC"]+=1
                        elif m.ALT[0] == "G":
                            basis["AT_CG"]
                    #figure out which of the 6 possible substs it was
            
    return basis

#%%
wholedict=dict()
for n in range(1,97):
    stringname="s1_" + str(n)
    try:
        thisfile=r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_VCF\s1_' + str(n) + r'_clean.VCF'
        wholedict[stringname]=count_fixed_mutations_s1(thisfile)
    except FileNotFoundError:
        continue
MA_s1_df=pd.DataFrame.from_dict(wholedict,orient="index")
MA_s1_df.to_csv(r'C:\Users\rmoge\GitHub\MuriSam\MA_s1_VCF\MA_s1_summary.csv') 
#%%
# get all sequence records for the specified genbank file
recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(recs))

# print annotations for each sequence record
for rec in recs:
	print(rec.annotations)
    
myrec=recs[0]

a=0;c=0;g=0;t=0
for nt in myrec:
    if nt == "A":
        a+=1
    elif nt == "C":
        c+=1
    elif nt == "G":
        g+=1
    elif nt == "T":
        t+=1
#a;c;g;t
'''
203606
67238
64720
207815
'''
syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(syn1recs))

# print annotations for each sequence record
for rec in syn1recs:
	print(rec.annotations)
    
syn1rec=syn1recs[0]
s1a=0;s1c=0;s1g=0;s1t=0
for nt in syn1rec:
    if nt == "A":
        s1a+=1
    elif nt == "C":
        s1c+=1
    elif nt == "G":
        s1g+=1
    elif nt == "T":
        s1t+=1
#s1a;s1c;s1g;s1t        
'''
405638
131103
127799
414269
'''

#%%
gb_file=r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb"
gb_record = SeqIO.read(open(gb_file,"r"), "genbank")
print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
print(repr(gb_record.seq))


#%%
print(myrec.features[2].location.start)
myrec[(myrec.features[2].location.start)]#this is how you access the 0th nt of the gb file...
#for the parallel stuff, we just need lengths
#%%
def gene_to_tag(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type=='gene':
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key")
                        answer[value]['index'].append(index)
                    else :
                        answer[value] = dict()
                        answer[value]['index'] = list()
                        answer[value]['index'].append(index)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['locus_tag']=feature.qualifiers['locus_tag']    
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
GD = gene_to_tag(myrec, "locus_tag", "gene")
#%%
# I think it may make more sense to build the dict all in one loop, rather than piping an updated dict from one function to the next, though. I attempt to do that in this chunk.
#We are gonna build a badass dict of dicts. Each locus tag will be a key. Its value will be a dict. In that dict, each key will be some kind of feature (e.g., length) and its value will be the corresponding value for THAT LOCUS.
def index_genbank_features(gb_record, feature_type, qualifier, qualifier_translation):
    max_gene_length = max([len(i) for i in gb_record.features[1:]])
    print(str(max_gene_length) + " is max gene length")
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else :
                        answer[value] = {}
                        answer[value]["index"] = index
            try:
                mygg=feature.qualifiers[qualifier_translation]#attempt to access the gene name for this locus tag
            except KeyError:
                mygg=feature.qualifiers[qualifier]#if there is no gene name, just call it by the locus tag
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene'] = mygg[0]
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene_length']=len(feature)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :        
                    answer[value]["gene_len_rel"]= (len(feature) / max_gene_length)
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
DD = index_genbank_features(myrec, "gene", "locus_tag", 'gene')
#%%
#Import mutation frequency data
def freq_list(infile=""):
    mydf = pd.read_csv(infile)
    return mydf          
            
syn3B_freqs = freq_list(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\mutation.frequencies_syn3B.csv")######THIS FILE INCLUDES SYNONYMOUS MUTATIONS, which means that your P-values may be overly conservative
syn3B_freqs_nosyn = freq_list(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\mutation.frequencies_no.syn_syn3B.csv")###in this file all the synonymous mutations have been deleted
#syn3B_freqs.set_index('freq',inplace=True)
syn3B_freqs_raw = tuple(syn3B_freqs_nosyn.loc[:, 'freq'].values)
#%%
roogeld = dict()
roogeld['JCVISYN3A_0001']= DD['JCVISYN3A_0001']; roogeld['JCVISYN3A_0002'] = DD['JCVISYN3A_0002']; roogeld['JCVISYN3A_0003'] = DD['JCVISYN3A_0003']

def simulate_mutations(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master
#%%
start = time.time()
smo = simulate_mutations((syn3B_freqs_raw),DD,reps = 100000)
end = time.time()
print(end - start)
#10000 simulations took 29.5 seconds. I could thus do 100 000 simulations in about 6 minutes, or 1 000 000 in 60 minutes.
#print(smo) 
###############################################################################################################
#Next steps:
#Compare simulated values to actual values. Get p values. Then, get Q values (FDR-adjusted P-values)
##could use the GENExPOP matrix to do that.
#%%
##Accessing the simulation list of values BY GENE NAME:
##homer = smo[GD['ftsZ']['locus_tag'][0]]
#o = np.count([i>=x for i in reps]), and I let reps = 10 000 in this simulation.
#for FtsZ, I get P < 0.0001      x = 4
#for FakA, I get P = 0.0016 x = 3    #note that these were calculated including synonymous mutations, which means the P-values may be overly conservative
    #tetM: P = 0.0019       x = 2.877
    #dnaN: P < 0.0001       x = 4
    #dnaA: P = 0.1388       x = 1
    #ptsG: P = 0.0264       x = 1.969
    #atpD: P < 0.0001        x = 3.052
    #rpoC: P = 0.0011       x = 4.056
    #rpoB: P = 0.3692      x = 1
    #amiF: P = 0.0023 x = 2.596 (synonymous pmsm NOT Included)
    #pyrG: P < 0.0001     x = 4
    
#recalculated without synonymous mutns in the simulation, this time with 100 000 reps
    
#%%  
    #these P-values were calculated without any synonymous mutations counted, 100 000 simulations
11*np.count_nonzero([i>=4 for i in smo[GD['ftsZ']['locus_tag'][0]]])######P = 0.00001 padj < 0.00011
11*np.count_nonzero([i>=4 for i in smo[GD['dnaN']['locus_tag'][0]]])######P < 0.00001 padj = 0.00022
11*np.count_nonzero([i>=4 for i in smo[GD['pyrG']['locus_tag'][0]]])#P < 0.00001      padj < 0.00011
11*np.count_nonzero([i>=1 for i in smo[GD['rpoB']['locus_tag'][0]]])#####P = 0.33659  padj = 1
11*np.count_nonzero([i>=4.056 for i in smo[GD['rpoC']['locus_tag'][0]]])##P = 0.00031 padj = 0.00484
11*np.count_nonzero([i>=0.456 for i in smo[GD['rpoD']['locus_tag'][0]]])##P = 0.16773 padj = 1
11*np.count_nonzero([i>=3 for i in smo[GD['fakA']['locus_tag'][0]]])######P = 0.00084 padj = 0.00748
11*np.count_nonzero([i>=1.969 for i in smo[GD['ptsG']['locus_tag'][0]]]) #P = 0.02367 padj = 0.25630 <-- So the Bonferroni corrected is no longer significant.
11*np.count_nonzero([i>=3.052 for i in smo[GD['atpD']['locus_tag'][0]]]) #P < 0.00001 padj < 0.00011
11*np.count_nonzero([i>=2.877 for i in smo[GD['tetM']['locus_tag'][0]]])##P = 0.00141 padj = 0.01463
11*np.count_nonzero([i>=2.596 for i in smo[GD['amiF']['locus_tag'][0]]])#P = 0.00121  padj = 0.01386

#TO GET BONFERRONI ADJUSTED P-VALUES, YOU MULTIPLIED THESE P-VALUES BY 11. You can also use the FDR and get Q-values, or use Dunn-Sidak, to get less conservative estimates.

#np.count_nonzero([i>=1 for i in smo[GD['dnaA']['locus_tag'][0]]])######P = 0.13024
#np.count_nonzero([i>=1.074 for i in smo[GD['amiF']['locus_tag'][0]]])#P = 0.04771


#############
#%%
#Now do 100 000 simulations for syn1.0
syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(syn1recs))

# print annotations for each sequence record
for rec in syn1recs:
	print(rec.annotations)
    
syn1rec=syn1recs[0]

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in syn1recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)

GD1 = gene_to_tag(syn1rec, "locus_tag", "gene")
DD1 = index_genbank_features(syn1rec, "gene", "locus_tag", 'gene')

syn1_freqs_nosyn = freq_list(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\mutation.frequencies_no.syn_syn1.0.csv")###in this file all the synonymous mutations have been deleted
#syn3B_freqs.set_index('freq',inplace=True)
syn1_freqs_raw = tuple(syn1_freqs_nosyn.loc[:, 'freq'].values)

#%%


start = time.time();print(start)
smo1 = simulate_mutations((syn1_freqs_raw),DD1,reps = 100000)
end = time.time()
print(end - start)

#%%  
    #these P-values were calculated without any synonymous mutations counted, 100 000 simulations
    #IGNORE ALL THE P VALUES---they are wrong. Only read the padj values. Divide them by 12 for the P value. Note that the padj values are from Bonferroni correction. To be less conservative, you could do Q value from FDR or Dunn-Sidak
12*np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]]])######P = 0.00001 padj < 0.00012
12*np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]]])######P < 0.00001 padj = 0.00108
12*np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]]])#P < 0.00001      padj < 0.02016
12*np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]]])#####P = 0.33659  padj = 0.94428
12*np.count_nonzero([i>=0.782 for i in smo1[GD1['rpoC']['locus_tag'][0]]])##P = 0.00031 padj = 1
12*np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]]])##P = 0.16773 padj = .72120
12*np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]]])######P = 0.00084 padj < 0.00012
12*np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]]]) #P = 0.02367 padj = 0.00336
12*np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]]]) #P < 0.00001 padj = 0.00552
12*np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]]])##P = 0.00141 padj = 0.08208
12*np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460']])##############padj = 0.33072
12*np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641']])#P = 0.00121  padj < 0.00012
#Dunn-Sidak
1-(1-(np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]]]))/100000)**12######P = 0.00001 padj < 0.00012
1-(1-(np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]]]))/100000)**12######P < 0.00001 padj = 001079465560347992
1-(1-(np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]]]))/100000)**12#P < 0.00001      padj < 0.019974760826477422
1-(1-(np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]]]))/100000)**12#####P = 0.33659  padj = 0.6260018785275774
1-(1-(np.count_nonzero([i>=0.782 for i in smo1[GD1['rpoC']['locus_tag'][0]]]))/100000)**12##P = 0.00031 padj = 0.6260018785275774
1-(1-(np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]]]))/100000)**12##P = 0.16773 padj = 0.5246868876744993
1-(1-(np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]]]))/100000)**12######P = 0.00084 padj < 0.00012
1-(1-(np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]]]))/100000)**12 #P = 0.02367 padj = 0.00335483042639817
1-(1-(np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]]]))/100000)**12 #P < 0.00001 padj = 0.005506055791773101
1-(1-(np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]]]))/100000)**12##P = 0.00141 padj = 0.07906148163292737
1-(1-(np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460']]))/100000)**12##############padj = 0.2849214088321784
1-(1-(np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641']]))/100000)**12#P = 0.00121  padj < 0.00012


#%%
#o = np.count([i>=x for i in reps])
#Let x = number of OBSERVED mutations for gene of interest.
#Let [reps] be the array of simulated mutation values for the gene of interest
#Then do this line of code for each gene (in order to do this, you need the GENExPOP matrix). Store a P-value in DD for each gene, where P = o / len(reps). We are asking, "What proportion of the time do we see AT LEAST as many mutations (as were actually observed) to occur BY CHANCE ALONE?"
#Next, you need to get Q values (FDR-adjusted P values)
##Rank the genes by P value
##Calculate Q as p_i * m / i, where i is the rank, p_i is the P value for the gene, and m is the total number of P values in your rankings

#That said--- the Q value and FDR are for GWAS where you are fishing for hits. If I just submit my candidate genes, such as ftsZ, I am NOT MAKING 500 TESTS. I am only testing the ~10 genes that I really think might be involved in adaptation.

#%%
#A next step could be to add the gene LENGTH to the dict. Later when youre worrying about dNdS, you would add the gene (DNA) SEQUENCE to the dict.OR MAYBE THE gff3 file would have an easier way to access teh DNA sdequence; worth looking at.
            
def looptest(dicter, reps=2):
    counter = 0
    while counter < reps:
        while True:
            w = random.choice(list(dicter)); r = random.random()
            print("w = " + str(w) + ", w len rel = " + str(dicter[w]['gene_len_rel']) + " and r = " + str(r))
            if dicter[w]['gene_len_rel'] >= r:
                print("sucksess")
                break
            else:
                continue
        counter += 1    
    
    
#looptest(DD)
    
    
    
#%%
def output_genes(gb_record, outfile=""):
    tempdf=pd.DataFrame()
    counter=0
    for (index,feature) in enumerate(gb_record.features):
        if feature.type == 'gene':
            if 'gene' in feature.qualifiers:
                for v in feature.qualifiers['gene']:
                    tempdf[str(feature.qualifiers['gene'][0])] = ""
            else:#~147 genes do not have assigned gene names and are known ONLY by the locus tag.
                for v in feature.qualifiers['locus_tag']:
                    tempdf[str(feature.qualifiers['locus_tag'][0])] = ""
            counter+=1
    print(counter)            

    tempdf.to_csv(index=True,path_or_buf=outfile)    
    return tempdf

gxp = output_genes(myrec,r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\gxp_headers_syn3B.csv")
#pp.pprint(gxp)
##
#%%
#Now output the simulation results
simulation_3B_100000 = pd.DataFrame.from_dict(smo)
simulation_1_100000 = pd.DataFrame.from_dict(smo1)

simulation_3B_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\simulations_syn3B_100000.csv")
simulation_1_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\simulations_syn1.0_100000.csv")
#%%
#You can import simulation results from a file.
start = time.time();print(start)
recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb", "genbank")]
myrec=recs[0]
syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]  
syn1rec=syn1recs[0]

GD = gene_to_tag(myrec, "locus_tag", "gene")
DD = index_genbank_features(myrec, "gene", "locus_tag", 'gene')
GD1 = gene_to_tag(syn1rec, "locus_tag", "gene")
DD1 = index_genbank_features(syn1rec, "gene", "locus_tag", 'gene')


smopd = pd.read_csv(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\simulations_syn3B_100000.csv")
smo1pd = pd.read_csv(r"C:\Users\rmoge\Box Sync\Mycoplasma\Strains\simulations_syn1.0_100000.csv")

smo = pd.DataFrame.to_dict(smopd)
smo1 = pd.DataFrame.to_dict(smo1pd)
end = time.time()
print(end - start)

#Now do the analyses on the imported simulation resulsts
12*np.count_nonzero([i>=4 for i in smo[GD['ftsZ']['locus_tag'][0]].values()])######P = 0.00001 padj < 0.00011
12*np.count_nonzero([i>=4 for i in smo[GD['dnaN']['locus_tag'][0]].values()])######P < 0.00001 padj = 0.00011
12*np.count_nonzero([i>=4 for i in smo[GD['pyrG']['locus_tag'][0]].values()])#P < 0.00001      padj = 0.00044
12*np.count_nonzero([i>=1 for i in smo[GD['rpoB']['locus_tag'][0]].values()])#####P = 0.33659  padj = 1
12*np.count_nonzero([i>=4.056 for i in smo[GD['rpoC']['locus_tag'][0]].values()])##P = 0.00031 padj = 0.00242
12*np.count_nonzero([i>=0.456 for i in smo[GD['rpoD']['locus_tag'][0]].values()])##P = 0.16773 padj = 1
12*np.count_nonzero([i>=3 for i in smo[GD['fakA']['locus_tag'][0]].values()])######P = 0.00084 padj = 0.00737
12*np.count_nonzero([i>=1.969 for i in smo[GD['ptsG']['locus_tag'][0]].values()]) #P = 0.02367 padj = 0.25532 <-- So the Bonferroni corrected is no longer significant.
12*np.count_nonzero([i>=3.052 for i in smo[GD['atpD']['locus_tag'][0]].values()]) #P < 0.00001 padj < 0.00011
12*np.count_nonzero([i>=2.877 for i in smo[GD['tetM']['locus_tag'][0]].values()])##P = 0.00141 padj = 0.01408
12*np.count_nonzero([i>=2.596 for i in smo[GD['amiF']['locus_tag'][0]].values()])#P = 0.00121  padj = 0.01463
12*np.count_nonzero([i>=2 for i in smo['JCVISYN3A_0430'].values()])##############padj = 0.00768

13*np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]].values()])######P = 0.00001 padj < 0.00012
13*np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]].values()])######P < 0.00001 padj = 0.00084
13*np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]].values()])#P < 0.00001      padj = 0.02112
13*np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]].values()])#####P = 0.33659  padj = 0.95400
13*np.count_nonzero([i>=0.782 for i in smo1[GD1['rpoC']['locus_tag'][0]].values()])##P = 0.00031 padj = 1
13*np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]].values()])##P = 0.16773 padj = .71580
13*np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]].values()])######P = 0.00084 padj < 0.00012
13*np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]].values()]) #P = 0.02367 padj = 0.00156
13*np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]].values()]) #P < 0.00001 padj = 0.00480
13*np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]].values()])##P = 0.00141 padj = 0.08544
13*np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460'].values()])##############padj = 0.32784
13*np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641'].values()])#P = 0.00121  padj < 0.00013
13*np.count_nonzero([i>=0.172 for i in smo1[GD1['rpsG']['locus_tag'][0]].values()])######P = 0.64831


#%%
def compare_recs(gb_record1,gb_record2):
    answer = dict()
    bads=list()
    count=0
    for (index, feature) in enumerate(gb_record1.features) :
        if feature.type=='CDS':
            if "locus_tag" in feature.qualifiers :
                lt1 = feature.qualifiers['locus_tag']
                num1 = lt1[0][-5:]
                #print(num1)
                try:
                    psq1 = feature.qualifiers['translation']
                except KeyError:
                    pass
                #try:
                 #   print(psq)
                #except UnboundLocalError:
                    #print("")
                for (index2,feature2) in enumerate(gb_record2.features):
                    if feature2.type=='CDS':
                        if "locus_tag" in feature2.qualifiers:
                            lt2 = feature2.qualifiers['locus_tag']
                            num2 = lt2[0][-5:]
                            if num1 == num2:
                                #print('match!')
                                count+=1
                                try:
                                    psq2 = feature2.qualifiers['translation']
                                    if psq1 != psq2:
                                        print('not the same')
                                        bads+=lt1
                                    else:
                                        print('same')
                                except:
                                    pass
                                       
                                       
    return bads

#test = index_genbank_features(myrec, "CDS", "gene")
cd=compare_recs(myrec, syn1rec)
print(cd)