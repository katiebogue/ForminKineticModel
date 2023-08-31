function [o]=polymerstats(formin,lt)
fh1_length=formin.fh1_length;
pp_index_vec=formin.pp_index_vec;
    [row,~] = find(lt.X1(:,1) == fh1_length & lt.X1(:,2) == pp_index_vec); %num of row corresponding to fH1 length AND location of PP
    o.p_occ1 = lt.X1(row,3); %all probabilities for that length at each PP location
    [row_0,~] = find(lt.X1(:,1) == fh1_length & lt.X1(:,2) == 1);
    o.p_occ1_0 = lt.X1(row_0,3)*ones([length(pp_index_vec) 1]); %probability at position 0
    o.p_r1 = lt.X1(row,4);
    
    [row,~] = find(lt.X2a(:,1) == fh1_length & lt.X2a(:,2) == pp_index_vec);
    o.p_occ2a = lt.X2a(row,3);
    [row_0,~] = find(lt.X2a(:,1) == fh1_length & lt.X2a(:,2) == 1);
    o.p_occ2a_0 = lt.X1(row_0,3)*ones([length(pp_index_vec) 1]);
    o.p_r2a = lt.X2a(row,4);
   
    [row,~] = find(lt.X2b(:,1) == fh1_length & lt.X2b(:,2) == pp_index_vec);
    o.p_occ2b = lt.X2b(row,3);
    [row_0,~] = find(lt.X2b(:,1) == fh1_length & lt.X2b(:,2) == 1);
    o.p_occ2b_0 = lt.X2b(row_0,3)*ones([length(pp_index_vec) 1]);
    o.p_r2b = lt.X2b(row,4);
    
    
     [row,~] = find(lt.X3a(:,1) == fh1_length & lt.X3a(:,2) == pp_index_vec);
     o.p_occ3a = lt.X3a(row,3);
     [row_0,~] = find(lt.X3a(:,1) == fh1_length & lt.X3a(:,2) == 1);
     o.p_occ3a_0 = lt.X3a(row_0,3)*ones([length(pp_index_vec) 1]);
     o.p_r3a = lt.X3a(row,4);
     
     [row,~] = find(lt.X3b(:,1) == fh1_length & lt.X3b(:,2) == pp_index_vec);
     o.p_occ3b = lt.X3b(row,3);
     [row_0,~] = find(lt.X3b(:,1) == fh1_length & lt.X3b(:,2) == 1);
     o.p_occ3b_0 = lt.X3b(row_0,3)*ones([length(pp_index_vec) 1]);
     o.p_r3b = lt.X3b(row,4);

end