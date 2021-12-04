#!/usr/bin/env python
# coding: utf-8

# In[1]:


def gathering(query):
    from bioservices import UniProt
    
    # Import BeautifulSoup, a package specialized for interpreting xml data
    from bs4 import BeautifulSoup

    # These packages need to be installed before being imported. If you have pip, then use
    # $ pip install beautifulsoup
    # $ pip install bioservices

    # Import itertools for efficient looping
    from itertools import groupby
    # Import numpy for efficient array/math functions
    import numpy as np
    from numpy import floor
    
    service = UniProt() 
    result_xml = service.search(query, frmt="xml")
    soup = BeautifulSoup(result_xml, 'html.parser') # xml "soup" object
    
    featureFH1 = soup.find_all('feature', description='FH1')
    
    # note the following code assumes there is one (and only one) annoted FH1 in this structure

    if len(featureFH1) == 0:
        print('No FH1 domain in this protein')



    beginPosition = int(featureFH1[0].find('location').find('begin').get('position'))
    endPosition = int(featureFH1[0].find('location').find('end').get('position'))

    lengthOfFH1 = endPosition-beginPosition+1
    
    if lengthOfFH1 == 0:
        print('error')
    #print(lengthOfFH1)
    
    soup_sequences = soup.find_all('sequence')
    soup_sequence = soup_sequences[-1].get_text();
    #print(soup_sequences)
    #print(soup_sequence)
     #script to display index and number of prolines in each PP sequence. Also returns 
        #vectors containing each information for ease of plotting


    soup_sequence = soup_sequence.replace("\n", "").replace("\r","") # gets rid of newline and line breaks in string 
    fh1_sequence = soup_sequence[beginPosition-1:endPosition] #specifiying FH1 domain
    fh1_sequence = fh1_sequence[::-1] #reverses string sequence
    #print(fh1_sequence)

    displayIndex = 0 # index used for poly_proline sequence
    index = 0 # regular indexing 

    pp_index_vec = [] # poly-proline vector for storing index
    pp_length_vec = [] # poly-proline vector for storing length

    fh1_length = len(fh1_sequence)
    fh1_length = float(fh1_length)
    #print(f'\nLength of entire sequence is {fh1_length}')

    seq = (groupby(fh1_sequence));


    for (k,g) in seq:
        print(k)
        length_seq = len(list(g)) # length of poly_proline sequence
        print(length_seq)
        if k=='P' and length_seq >3: # for indexing, refer to report
            if length_seq%2 == 0: 
                displayIndex = floor ((length_seq) / 2 ) - 1 + index 
                pp_index_vec.append(displayIndex)
                pp_length_vec.append(length_seq)
            else:
                displayIndex = floor ((length_seq) / 2 ) + index
                pp_index_vec.append(displayIndex)
                pp_length_vec.append(length_seq)
        index += length_seq
        
    pp_length_vec = [float(i) for i in pp_length_vec]
        
    #print(f'\nPoly_proline index vector:{pp_index_vec}')
    #print(f'\nPoly_proline length vector: {pp_length_vec}')
    
    return fh1_length, pp_index_vec, pp_length_vec


# In[3]:


#test case
# Q24120 = Capu
gathering('Q24120')


# In[ ]:





# In[ ]:




