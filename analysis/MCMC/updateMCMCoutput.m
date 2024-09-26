function updateMCMCoutput(matfileloc,mem)
%UPDATEMCMCOUTPUT modified variables in target MCMCParamfit output (.mat
%file) to remove empty elements arising from preallocation as well as
%convert things into non log scale if needed.
%
%   UPDATEMCMCOUTPUT(matfileloc) modify variables in matfileloc, stepping
%   through arrays assuming a memory of 15.5GB
%
%   UPDATEMCMCOUTPUT(matfileloc, mem) modify variables in matfileloc, stepping
%   through arrays based on specified mem
%
%   Inputs:
%       matfileloc : path to the target .mat file with MCMC output
%       mem        : GB of RAM to guide how much of each array is loaded
%           into memory at once (default is 15.5)
%       
%   Loads file as a matfile and then does the following:
%       - removes trailing zeros in paccept_matrix
%       - removes trailing zeros in ksvals
%       - removes trailing zeros in paramHistCounts_matrices_nts
%       - removes empty sheets in paramHistCounts_matrices
%       - creates logparams_all_trun, which is logparams_all without
%       traliling zeros
%       - creates corruptindex, which contains the locations of any
%       corruptions in logparams_all, for which logparams_all_trun contains
%       NaN
%       - if it does not exist, creates parameters_all, with values from
%       logparams_all_trun not in log space; or, truncates parameters_all
%       - creates containsInf, which is a series of bools saying whether
%       there is a +/- Inf value in the parameters
%       - creates updatedone and sets it to true
%       - sets logparams_all to an empty array
% 
% See also VISUALIZEPOSTERIORS, MATFILE, MCMCPARAMFIT.
arguments
    matfileloc
    mem=15.5 % GB RAM
end

    m = matfile(matfileloc,'Writable',true);

    if ~isempty(who(m,'updatedone'))
        if m.updatedone
            disp("update already done, not updating")
            return
        end
    end

    if ~isempty(who(m,'paccept_matrix'))
        m.paccept_matrix( ~any(m.paccept_matrix,2),:) = [];
        disp("shortened paccept_matrix")
    end
    if ~isempty(who(m,'ksvals'))
        m.ksvals( ~any(m.ksvals,2),:) = [];
        disp("shortened ksvals")
    end
    if ~isempty(who(m,'paramHistCounts_matrices'))
        m.paramHistCounts_matrices_nts( ~any(m.paramHistCounts_matrices_nts,2)) = [];
        disp("shortened paramHistCounts_matrices_nts")
        temp=m.paramHistCounts_matrices(:,:,1:length(m.paramHistCounts_matrices_nts));
        m.paramHistCounts_matrices=temp;
        disp("shortened paramHistCounts_matrices")
        clear temp
    end

    [nrows,ncols]=size(m,'logparams_all');
    
    maxrow=m.nt;
    if all(m.logparams_all(maxrow,:)==0)
       maxrow=maxrow-1;
       while maxrow>0 && all(m.logparams_all(maxrow,:)==0)
           maxrow=maxrow-1;
       end
    else
       maxrow=maxrow+1;
       len=size(m,'logparams_all',1);
       while maxrow<=len && ~all(m.logparams_all(maxrow,:)==0)
           maxrow=maxrow+1;
       end
       maxrow=maxrow-1;
    end
    disp("found maxrow")

    m.logparams_all_trun(maxrow,ncols)=0;

    corruptindex=[];
    i=1;

    STEP=uint64(mem*(1024^3)/(ncols*8));
    disp("starting truncation")
    while i+STEP<=maxrow
        try
            x=m.logparams_all(i:i+STEP,:);
            m.logparams_all_trun(i:i+STEP,:)=x;
        catch
            disp("unable to load so trying to loop through each")
            for j=i:uint64(i+STEP)
                try
                    x=m.logparams_all(j,:);
                    m.logparams_all_trun(j,:)=x;
                catch
                    m.logparams_all_trun(j,:)=NaN;
                    corruptindex=[corruptindex; j];
                end
            end
        end
        i=i+STEP+1;
        disp(i)
    end
    disp("finished initial loop")
    ind=i:maxrow;
    try
        x=m.logparams_all(ind,:);
        m.logparams_all_trun(ind,:)=x;
    catch
        for j=ind
            
            try
                x=m.logparams_all(j,:);
                m.logparams_all_trun(j,:)=x;
            catch
                m.logparams_all_trun(j,:)=NaN;
                corruptindex=[corruptindex; j];
            end
        end
    end
    m.corruptindex=corruptindex;

    disp("logparams_all_trun successfully made")

    if isempty(who(m,'parameters_all'))||length(m.parameters_all)==1
        memarraysize=16*(1024^3)/(2*8); %limit to the largest possible in memory array size for 2 rows (so hist3 works)
        thirdarraysize=ceil(maxrow/3)*2;
        if memarraysize>thirdarraysize
            maxarraysize=thirdarraysize;
        else
            maxarraysize=memarraysize;
        end
        if maxarraysize<STEP
            m.parameters_all=m.logparams_all_trun(maxrow-maxarraysize+1:maxrow,:);
            nparams=m.nparams;
            m.containsInf=zeros(nparams,1);
            len=size(m,'parameters_all',1);
            for i=1:nparams
                if max(m.parameters_all(1:len,i))==Inf || min(m.parameters_all(1:len,i))==-Inf
                    m.containsInf(i,1)=1;
                end
            end
            disp("made containsInf")
        else
            disp("unable to do 1/3, so trying based on size")
            nparams=m.nparams;
            m.containsInf=zeros(nparams,1);
            m.parameters_all(maxarraysize,ncols)=0;
            minrow=uint64(maxrow-maxarraysize+1);
            %disp("starting loops to save parameters_all")
            i=1;
            while i+STEP<=maxarraysize
                disp(i)
                lowerbound=uint64(minrow+i-1);
                %disp(lowerbound)
                upperbound=uint64(minrow+STEP+i-1);
                %disp(upperbound)
                x=m.logparams_all_trun(lowerbound:upperbound,:);
                disp("saved x")
                m.parameters_all(i:uint64(i+STEP),:)=x;
                for j=1:nparams
                    if m.containsInf(j,1) || max(x(:,j))==Inf || min(x(:,j))==-Inf
                        m.containsInf(j,1)=1;
                    end
                    disp("made containsInf")
                end
                i=i+STEP+1;
            end
            x=m.logparams_all_trun(minrow+i-1:maxrow,:);
            m.parameters_all(i:maxarraysize,:)=x;
            for j=1:nparams
                if m.containsInf(j,1) || max(x(:,j))==Inf || min(x(:,j))==-Inf
                    m.containsInf(j,1)=1;
                end
                %disp("made containsInf")
            end
            clear x
        end
    else
        nparams=m.nparams;
        m.containsInf=zeros(nparams,1);
        len=size(m,'parameters_all',1);
        for i=1:nparams
            if max(m.parameters_all(1:len,i))==Inf || min(m.parameters_all(1:len,i))==-Inf
                m.containsInf(i,1)=1;
            end
        end
        disp("made containsInf")
    end
    disp("made parameters_all")

    m.logparams_all=[];
    disp("removed logparams_all")

    m.updatedone=1;
end