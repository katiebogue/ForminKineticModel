function [pocc_dict, pr_dict,pocc_ext,pr_ext]=readlookup(nFil,iFil,file,maxL)
% READLOOKUP  Reads in a lookup table .txt file and identifies occlusion
% probabilities and pr vectors. 
    %
    %   [pocc_dict, pr_dict,pocc_ext,pr_ext]=
    %   READLOOKUP(nFil,iFil,file,maxL) generates dictionaries for pocc and
    %   pr for the iFil-th filament of nFil filaments with a maximum length
    %   of maxL by reading in the lookup table in file
    %
    %   Inputs:
    %         nFil : total number of filaments (1 for single, 2 for dimer
    %                and double)
    %         iFil : filament number 0 or 1 (for dimer and double)
    %         file : lookup table file
    %         maxL : max filament length before extrapolation
    %   
    %   Outputs:
    %         pocc_dict : dictionary with filament length as keys and
    %                     arrays containing occlusion probabilities for
    %                     each iSite (in order) as values
    %         pr_dict   : dictionary with filament length as keys and
    %                     arrays containing pr at the barbed end for each
    %                     iSite (in order) as values
    %         pocc_ext  : scatteredInterpolant of pocc values (access via
    %                     pocc_ext(FH1length, iSiteLoc)
    %         pr_ext    : scatteredInterpolant of pr values (access via
    %                     pocc_ext(FH1length, iSiteLoc)
    %
    %   See also LOOKUPTABLE.
    
    numErr=0;
    tab = dlmread(file);
    pocc_dict = dictionary;
    pr_dict = dictionary;
    length_vec = [];
    iSite_vec = [];
    p_occ_all = [];
    p_r_all = [];
    
    
    for iFilLen=1:maxL
        siteCount=nFil*iFilLen;
        poccs = [];
        prs = [];
        for iSite =1:iFilLen
            length_vec = [length_vec iFilLen]; % list of filament lengths
            iSite_vec=[iSite_vec iSite]; % all iSites used
            p_occ = tab(iFilLen, 16 + 2*(siteCount +1) + 7*(iSite - 1) + iFil*(6 + 9*iFilLen + 2 + nFil + nFil));
            p_r = tab(iFilLen, 19 + 2*(siteCount +1) + 7*(iSite - 1) + iFil*(6 + 9*iFilLen + 2 + nFil + nFil));
            if p_occ > 1 || p_r >1
                numErr = numErr+1;
            end
            poccs=[poccs,p_occ];
            prs=[prs,prs];
            p_occ_all=[p_occ_all p_occ];
            p_r_all=[p_r_all p_r];
        end
        pocc_dict{iFilLen}=poccs;
        pr_dict{iFilLen}=prs;
    end
    
    X = [length_vec; iSite_vec; p_occ_all; p_r_all]' ;
    pocc_ext = scatteredInterpolant(X(:,1),X(:,2),X(:,3),'linear','nearest');
    pr_ext = scatteredInterpolant(X(:,1),X(:,2),X(:,4),'linear','nearest');
    
end