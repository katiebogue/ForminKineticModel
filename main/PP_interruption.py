#!/usr/bin/env python
# coding: utf-8

# In[4]:



# coding: utf-8

# In[35]:


def gathering_int(query,min_PP_length):
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
    from numpy import ceil
    
    service = UniProt() 
    result_xml = service.search(query, frmt="xml")
    soup = BeautifulSoup(result_xml, 'html.parser') # xml "soup" object
    
    featureFH1 = soup.find_all('feature', description='FH1')
    
    # note the following code assumes there is one (and only one) annoted FH1 in this structure!

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
    pp_index_vec2 = [] # poly-proline vector for storing index
    pp_length_vec2 = [] 
    
    fh1_length = len(fh1_sequence)
    fh1_length = float(fh1_length)
    #print(f'\nLength of entire sequence is {fh1_length}')


    seq = (groupby(fh1_sequence)); # group string by letter
    length_seq=[]
    letter_seq=[]
    seq_index = 0

    
    for (k,g) in seq:
        length_seq.append(len(list(g)))
        letter_seq.append(k) 
                      
    #print(length_seq, letter_seq)
    #print(len(letter_seq))

    #for loop to traverse each letter of the sequence. Need length-2 because last letter can't search for two ahead of it
    #if the letter is P, and the next length after that is 1 and not P, and the letter two ahead of it is P
    #find the total length of that sequence
    #if the total length of that sequence is even, pick the P closest to C-terminal as the "middle"
    #if the total length of that sequence is odd, pick the P in the middle
    #return the length of fh1, the index of the PP sequence, and the length of the PP sequence
    #now we also have to account for all the P's without interruptions but it'll count the ones with interruptions twice (PPPPXP)
    #if the letter is P and the length is greater than 4 and the letter two before or after it isn't P
    #find the total length of that sequence
    #if the total length of that sequence is even, pick the P closest to C-terminal as the "middle"
    #if the total length of that sequence is odd, pick the P in the middle
    #return the length of fh1, the index of the PP sequence, and the length of the PP sequence

    n = 0
    while n < len(letter_seq):
    #for n in range (0, len(letter_seq)): 
    #for n in range (0,len(letter_seq)):
    #print(letter_seq[n] , length_seq[n])
    #print(n)
        if n < len(letter_seq)-2 and letter_seq[n] == 'P' and length_seq[n+1] ==1 and letter_seq[n+2] == 'P':
            total_length_PRM = length_seq[n] +  length_seq[n+1] + length_seq[n+2]
            #print('index', n, 'total length', total_length_PRM, 'length of sequence', length_seq[n])
            if total_length_PRM>= min_PP_length:
                displayIndex = ceil ((total_length_PRM) / 2 )  + seq_index 
                pp_index_vec.append(displayIndex)
                pp_length_vec.append(total_length_PRM)
                n = n+2;
                seq_index += total_length_PRM
            else:
                seq_index += length_seq[n]
                
                
    #4P without interruption
        elif letter_seq[n] == 'P' and length_seq[n] > 5: 
        #and letter_seq[n-2] != 'P' and letter_seq[n+2] !='P':
        #problem is if the PP is at the end of the beginning. can't look two behind or two ahead. might give an error?
        #total_length2 = length_seq[n] 
        #print(n, total_length2, length_seq[n])
            #if length_seq[n]%2 == 0: 
                displayIndex = ceil ((length_seq[n]) / 2 ) + seq_index
                pp_index_vec.append(displayIndex)
                pp_index_vec2.append(displayIndex)
                pp_length_vec.append(length_seq[n])
                pp_length_vec2.append(length_seq[n])
                seq_index +=length_seq[n]
        else:
            seq_index += length_seq[n]
        #print(n, length_seq[n]) 
    
        n += 1
        
        
     
    pp_length_vec = [float(i) for i in pp_length_vec]
    pp_length_vec2 = [float(i) for i in pp_length_vec2]
    
    return fh1_length, pp_index_vec, pp_length_vec, #pp_index_vec2, pp_length_vec2

    #return fh1_length, pp_index_vec2, pp_length_vec2, #pp_index_vec2, pp_length_vec2

    #return 24,24,24,24,24


# In[36]:


#gathering('Q24120')


# In[38]:


#gathering('Q9Y4D1')


# In[ ]:




