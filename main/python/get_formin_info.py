#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 13:38:32 2023

@author: Katiebogue

This is a script that reads a formin sequence (either from text or from UniProt) 
and determines the FH1 domain and PRMs.

arguments for main:
    1) UniProt ID ("-1" if using a sequence directly)
    2) Sequence NT to CT ("-1" if using UniProt)
    3) Minimum PRM length (including interruptions)
    4) Number of allowed interruptions 
        (counts once for each amino acid-- i.e. if the max int len is 2, AA is  acceptable but counts as 2 interruptions
    5) Maximum interruption length
    6) Minimum number of Ps 
    7) FH1 NT definition options:
        1 - first instance of PRM with at least 4 Ps with max 1 interruption of length 1
        2 - first instance of PRM (as defined by args 3 and 4)
        3 - Uniprot defined FH1 start (or option 1 if no FH1)
        4 - Uniprot defined FH1 start (or option 2 if no FH1)
        5 - Start of sequence (for input sequences)
    8) FH1 CT definition options:
        1 - Uniprot defined FH2 start
        2 - Uniprot defined FH1 end (or option 1 if no FH1)
        3 - End of sequence (for input sequences)
        (opt 3 is automatically selected if a sequence is provided in arg2)  
    9) PRM scan options:
        1 - Search for PRMs starting from the FH2 (CT) (default)
        2 - Search for PRMs starting from the FH1 NT

These packages need to be installed. If you have pip, then use
     $ pip install beautifulsoup
     $ pip install bioservices
"""
from bioservices import UniProt # to access UniProt 
from bs4 import BeautifulSoup # package specialized for interpreting xml data
from itertools import groupby # for efficient looping
import numpy as np # efficient array/math functions
from numpy import ceil # efficient array/math functions

#from sys import argv # allow argument input


def get_uniprot_info(uniprotID):
    service = UniProt()
    result_xml = service.search(uniprotID, frmt="xml")
    soup = BeautifulSoup(result_xml, 'html.parser') # xml "soup" object
    
    soup_sequence = soup.find_all('sequence') #finding full formin amino acid sequence
    soup_sequence = soup_sequence[-1].get_text();
    
    featureFH1 = soup.find_all('feature', description='FH1')
    featureFH2 = soup.find_all('feature', description='FH2') #adding FH2 domain to define end of FH1 domain

    if len(featureFH1)==0:
        FH1_start=-1
        FH1_end=-1
    else:
        FH1_start = int(featureFH1[0].find('location').find('begin').get('position'))
        FH1_end = int(featureFH1[0].find('location').find('end').get('position'))
    FH2_start = int(featureFH2[0].find('location').find('begin').get('position')) #end at FH2 domain

    return soup_sequence.upper(), FH1_start,FH1_end,FH2_start

def get_first_PRM(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP):
    seq=seq.upper()
    i=-1
    loc=-1
    nP=0
    lenPRM=0
    nInt=0
    lenInt=0
    while (i+1)< len(seq):
        i+=1
        #print(seq[i])
        if seq[i]=='P':
            lenInt=0
            if loc==-1:
                loc=float(i)
                nInt=0
            nP+=1
            lenPRM+=1
            if(nP>=min_nP and lenPRM>=min_PP_length):
                #print(nP)
                #print(lenPRM)
                return loc
        else:
            if loc ==-1:
                continue
            nInt+=1
            lenInt+=1
            if(nInt>nInterruptions or lenInt>lenInterruptions):
                loc=-1
                nInt=0
                lenInt=0
                lenPRM=0
                nP=0
            else:
                lenPRM+=1
    return loc            

def get_FH1_start(seq,NT_opt,min_PP_length,nInterruptions,lenInterruptions,min_nP,UP_FH1_start):
    if NT_opt==1:
        # FH1 starts at first PRM w/ min length 4, max 1 int of length 1, at least 4 Ps
        FH1_start=get_first_PRM(seq,4,1,1,4) 
    elif NT_opt==2:
        FH1_start=get_first_PRM(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP)
    elif NT_opt==3:
        if UP_FH1_start!=-1:
            FH1_start==UP_FH1_start
        else:
            FH1_start=get_first_PRM(seq,4,1,1,4) 
    elif NT_opt==4:
        if UP_FH1_start!=-1:
            FH1_start==UP_FH1_start
        else:
            FH1_start=get_first_PRM(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP)
    elif NT_opt==5:
        FH1_start=1
    
    return FH1_start

def get_FH1_end(seq,CT_opt,UP_FH1_end,UP_FH2_start):
    if CT_opt==1:
        FH1_end=UP_FH2_start
    elif CT_opt==2:
        if UP_FH1_end != -1:
            FH1_end=UP_FH1_end+1
        else:
            FH1_end=UP_FH2_start
    elif CT_opt==3:
        FH1_end=len(seq)
    
    return FH1_end

def get_FH1_seq(seq,NT_opt,CT_opt,\
                min_PP_length,nInterruptions,lenInterruptions,min_nP,\
                UP_FH1_start,UP_FH1_end,UP_FH2_start):
    
    FH1_start=get_FH1_start(seq,NT_opt,min_PP_length,nInterruptions,lenInterruptions,min_nP,UP_FH1_start)
    FH1_end=get_FH1_end(seq,CT_opt,UP_FH1_end,UP_FH2_start)
    FH1_sequence=seq[int(FH1_start-1):int(FH1_end)]
    
    return FH1_sequence
    
def get_Ps(seq):
    P_index_vec=[]
    seq=seq.upper()
    
    i=0
    while i<len(seq):
        if seq[i]=='P':
            P_index_vec.append(float(i))
        i+=1
    
    return P_index_vec
    
def get_PRMs(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP):
    seq=seq.upper()
    
    pp_index_vec = [] # poly-proline vector for storing index
    pp_length_vec = [] # poly-proline vector for storing length
    pp_index_vec_start = []
    P_index_vec=[]

    i=-1
    loc=-1
    nP=0
    nInt=0
    lenInt=0
    lastP=0
    p2=0 # next P after any interruptions
    while (i+1)< len(seq):
        i+=1
        #print(i)
        #print(seq[i])
        if seq[i]=='P':
            if (float(i)+1) not in P_index_vec:
                P_index_vec.append(float(i)+1)
            lenInt=0
            if loc==-1:
                loc=float(i)
                nInt=0
                p2=0
            if (nInt>=1 and p2==0):
                #print('assigned p2')
                # if there has been an int and p2 has not been assigned
                p2=int(i)
            nP+=1
            lastP=float(i)
            if i+1==len(seq):
                lenPRM=lastP-loc+1
                if(nP>=min_nP and lenPRM>=min_PP_length):
                    pp_index_vec_start.append(1+loc)
                    pp_length_vec.append(lenPRM)
                    #print(lenPRM)
                    PRM_loc=ceil(lenPRM/2) + loc
                    #print(PRM_loc)
                    pp_index_vec.append(PRM_loc)
        else:
            if loc ==-1:
                continue
            nInt+=1
            lenInt+=1
            if(nInt>nInterruptions or lenInt>lenInterruptions):
                lenPRM=lastP-loc+1
                if(nP>=min_nP and lenPRM>=min_PP_length):
                    pp_index_vec_start.append(1+loc)
                    pp_length_vec.append(lenPRM)
                    #print(lenPRM)
                    PRM_loc=ceil(lenPRM/2) + loc
                    #print(PRM_loc)
                    pp_index_vec.append(PRM_loc)
                elif p2!=0:
                    i=p2-1
                loc=-1
                nInt=0
                lenInt=0
                nP=0
            
    return pp_index_vec, pp_length_vec, pp_index_vec_start, P_index_vec
    

def main(a1,a2,a3,a4,a5,a6,a7,a8,a9):
    
    #input arguments as global variables:
    UniProtID=str(a1)
    Seq=str(a2) # Input sequence
    min_PP_length=float(a3)
    nInterruptions=float(a4)
    lenInterruptions=float(a5)
    min_nP=float(a6)
    FH1_NT_opt=float(a7)
    FH1_CT_opt=float(a8)
    PRM_search_opt=float(a9)

    # Get the full sequence
    if UniProtID != "-1":
        full_sequence,UP_FH1_start,UP_FH1_end,UP_FH2_start = get_uniprot_info(UniProtID)
    elif Seq != "-1":
        full_sequence=Seq.upper().replace("\n", "").replace("\r","")
        # gets rid of newline and line breaks in string, makes string uppercase
        FH1_CT_opt=3 #always use end of sequence
        UP_FH1_start=-1
        UP_FH1_end=-1
        UP_FH2_start=-1
    else:
        raise Exception("No UniProtID or sequence was provided")
    full_sequence=full_sequence.upper()
    #print(full_sequence)
    FH1_sequence = get_FH1_seq(full_sequence,FH1_NT_opt,FH1_CT_opt,\
                 min_PP_length,nInterruptions,lenInterruptions,min_nP,\
                 UP_FH1_start,UP_FH1_end,UP_FH2_start)
    #print(FH1_sequence)
    FH1_sequence = FH1_sequence.upper()    
    fh1_length=len(FH1_sequence)
    
    if PRM_search_opt==2:
        seq=FH1_sequence
        pp_index_vecr, pp_length_vec, pp_index_vec_startr, P_index_vecr = get_PRMs(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP)
        pp_index_vec=[]
        pp_index_vec_start=[]
        P_index_vec=[]
        for index in range(len(pp_index_vecr)):
            # put back in reference to dist to FH2
            if (pp_length_vec[index]%2) == 0:
                pp_index_vec.append(fh1_length-pp_index_vecr[index])
            else:
                pp_index_vec.append(fh1_length-pp_index_vecr[index]+1)
            pp_index_vec_start.append(fh1_length-(pp_index_vec_startr[index]+pp_length_vec[index]-1)+1)
        for x in P_index_vecr:
            P_index_vec.append(fh1_length-x+1)
        # reverse lists
        pp_index_vec=pp_index_vec[::-1]
        pp_length_vec=pp_length_vec[::-1]
        pp_index_vec_start=pp_index_vec_start[::-1]
        P_index_vec=P_index_vec[::-1]

    else:
        seq=FH1_sequence[::-1] # reverse sequence
        pp_index_vec, pp_length_vec, pp_index_vec_start, P_index_vec = get_PRMs(seq,min_PP_length,nInterruptions,lenInterruptions,min_nP)

    
        
    
    #print(fh1_length,'\n', pp_index_vec, '\n', pp_length_vec, '\n', P_index_vec,'\n', pp_index_vec_start)
    return fh1_length, pp_index_vec, pp_length_vec, P_index_vec, pp_index_vec_start 

#main()