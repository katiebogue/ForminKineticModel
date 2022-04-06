
# coding: utf-8

# In[35]:


# def gathering(query,min_PP_length):
#     from bioservices import UniProt
    
#     # Import BeautifulSoup, a package specialized for interpreting xml data
#     from bs4 import BeautifulSoup

#     # These packages need to be installed before being imported. If you have pip, then use
#     # $ pip install beautifulsoup
#     # $ pip install bioservices

#     # Import itertools for efficient looping
#     from itertools import groupby
#     # Import numpy for efficient array/math functions
#     import numpy as np
#     from numpy import floor  
    
#     service = UniProt() 
#     result_xml = service.search(query, frmt="xml")
#     soup = BeautifulSoup(result_xml, 'html.parser') # xml "soup" object
    
#     featureFH1 = soup.find_all('feature', description='FH1')
    
#     # note the following code assumes there is one (and only one) annoted FH1 in this structure!

#     if len(featureFH1) == 0:
#         print('No FH1 domain in this protein')



#     beginPosition = int(featureFH1[0].find('location').find('begin').get('position'))
#     endPosition = int(featureFH1[0].find('location').find('end').get('position'))

#     lengthOfFH1 = endPosition-beginPosition+1
    
#     if lengthOfFH1 == 0:
#         print('error')
#     #print(lengthOfFH1)
    
#     soup_sequences = soup.find_all('sequence')
#     soup_sequence = soup_sequences[-1].get_text();
#     #print(soup_sequences)
#     #print(soup_sequence)
#      #script to display index and number of prolines in each PP sequence. Also returns 
#         #vectors containing each information for ease of plotting


#     soup_sequence = soup_sequence.replace("\n", "").replace("\r","") # gets rid of newline and line breaks in string 
#     fh1_sequence = soup_sequence[beginPosition-1:endPosition] #specifiying FH1 domain
#     fh1_sequence = fh1_sequence[::-1] #reverses string sequence
#     #print(fh1_sequence)

#     displayIndex = 0 # index used for poly_proline sequence
#     index = 0 # regular indexing 

#     pp_index_vec = [] # poly-proline vector for storing index
#     pp_length_vec = [] # poly-proline vector for storing length

#     fh1_length = len(fh1_sequence)
#     fh1_length = float(fh1_length)
#     #print(f'\nLength of entire sequence is {fh1_length}')


#     seq = (groupby(fh1_sequence)); # group string by letter

#     for (k,g) in seq:
#         length_seq = len(list(g)) # length of poly_proline sequence
#         if k=='P' and length_seq > min_PP_length: # for indexing, refer to report
#             if length_seq%2 == 0: 
#                 displayIndex = floor ((length_seq) / 2 ) - 1 + index 
#                 pp_index_vec.append(displayIndex)
#                 pp_length_vec.append(length_seq)
#             else:
#                 displayIndex = floor ((length_seq) / 2 ) + index
#                 pp_index_vec.append(displayIndex)
#                 pp_length_vec.append(length_seq)
#         index += length_seq
        
#     pp_length_vec = [float(i) for i in pp_length_vec]
        
#     #print(f'\nPoly_proline index vector:{pp_index_vec}')
#     #print(f'\nPoly_proline length vector: {pp_length_vec}')
    
#     return fh1_length, pp_index_vec, pp_length_vec

#%%

def gathering(query,min_PP_length,standard):


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
    featureFH2 = soup.find_all('feature', description='FH2') #adding FH2 domain to define end of FH1 domain
    
    # note the following code assumes there is one (and only one) annoted FH1 in this structure!
    
    soup_sequences = soup.find_all('sequence') #finding full formin amino acid sequence
    soup_sequencee = soup_sequences[-1].get_text();
    soup_sequence=[]
    j=0
    while j <= len(soup_sequencee)-1:
        if soup_sequencee[j] == 'P':
            soup_sequence.append("P")
            j+=1
        else:
            soup_sequence.append("X")
            j+=1
    #print(soup_sequence)
    soup_sequence= ''.join(soup_sequence)
    #print(soup_sequence)

    find_1= soup_sequence.find('PXPPP')
    find_2= soup_sequence.find('PPXPP')
    find_3= soup_sequence.find('PPPXP')
    find_4= soup_sequence.find('PPPP')
    
    find_all= [find_1,find_2,find_3,find_4]
    #print(find_all)
    for i in range(0,3):
        if find_all[i]>0:
            find_all[i]=find_all[i]
        else:
            find_all[i]= len(soup_sequence)+5
    # note the following code assumes there is one (and only one) annoted FH1 in this structure!

    if standard == 'Y': #set in wholepkg
        beginPosition = min(find_all) #start at first P of a series of at least four Ps with a max of 1 interruption
        endPosition = int(featureFH2[0].find('location').find('begin').get('position')) #end at FH2 domain
    elif len(featureFH1) == 0: #if uniprot does not have fh1 domain already defined
        beginPosition = min(find_all) #start at first P of a series of at least four Ps with a max of 1 interruption
        endPosition = int(featureFH2[0].find('location').find('begin').get('position')) #end at FH2 domain
    else: #if uniprot has fh1 domain already defined
        beginPosition = int(featureFH1[0].find('location').find('begin').get('position'))
        endPosition = int(featureFH1[0].find('location').find('end').get('position'))

    #print(beginPosition, endPosition)
    
    lengthOfFH1 = endPosition-beginPosition+1
    
    if lengthOfFH1 == 0:
        print('error')
    #print(lengthOfFH1)
    
    
    # print(soup_sequences)
    # print("  ")
    # print(soup_sequence)
    # print("  ")

     #script to display index and number of prolines in each PP sequence. Also returns 
        #vectors containing each information for ease of plotting
    
    
    soup_sequence = soup_sequence.replace("\n", "").replace("\r","") # gets rid of newline and line breaks in string 
    fh1_sequence = soup_sequence[beginPosition-1:endPosition] #specifiying FH1 domain
    fh1_sequence = fh1_sequence[::-1] #reverses string sequence
   # print(fh1_sequence)
   
    P_index_vec=[]
    j=0
    while j <= len(fh1_sequence)-1:
        if fh1_sequence[j] == 'P':
            P_index_vec.append(float(j))
            j+=1
        else:
            j+=1
   
    displayIndex = 0 # index used for poly_proline sequence
    index = 0 # regular indexing 
    
    pp_index_vec = [] # poly-proline vector for storing index
    pp_length_vec = [] # poly-proline vector for storing length
    pp_index_vec_start = []
    
    fh1_length = len(fh1_sequence)
    fh1_length = float(fh1_length)
    #print(f'\nLength of entire sequence is {fh1_length}')
    
    
    seq = (groupby(fh1_sequence)); # group string by letter
    
    for (k,g) in seq:
        length_seq = len(list(g)) # length of poly_proline sequence
        if k=='P' and length_seq >= min_PP_length: # for indexing, refer to report
            if length_seq%2 == 0: 
                displayIndex = floor ((length_seq) / 2 ) - 1 + index 
                pp_index_vec.append(displayIndex)
                pp_length_vec.append(length_seq)
                pp_index_vec_start.append(float(index))
            else:
                displayIndex = floor ((length_seq) / 2 ) + index
                pp_index_vec.append(displayIndex)
                pp_length_vec.append(length_seq)
                pp_index_vec_start.append(float(index))
        index += length_seq
        
    pp_length_vec = [float(i) for i in pp_length_vec]
        
    #print(f'\nPoly_proline index vector:{pp_index_vec}')
    #print(f'\nPoly_proline length vector: {pp_length_vec}')
    
    return fh1_length, pp_index_vec, pp_length_vec, P_index_vec, pp_index_vec_start 

# In[36]:


gathering('Q24120',3, 'Y')
#%%
#gathering('Q9VSQ0',4)



# In[38]:


#gathering('Q9Y4D1')

