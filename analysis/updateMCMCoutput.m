function updateMCMCoutput(matfileloc,mem)
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


    [nrows,ncols]=size(m,'params_all');

    maxrow=m.nt;
    if all(m.params_all(maxrow,:)==0)
       maxrow=maxrow-1;
       while maxrow>0 && all(m.params_all(maxrow,:)==0)
           maxrow=maxrow-1;
       end
    else
       maxrow=maxrow+1;
       len=size(m,'params_all',1);
       while maxrow<=len && ~all(m.params_all(maxrow,:)==0)
           maxrow=maxrow+1;
       end
       maxrow=maxrow-1;
    end
    disp("found maxrow")

    m.params_all_trun(maxrow,ncols)=0;
    corruptindex=[];
    i=1;

    newlog=0;
    if isempty(who(m,'logparams_all'))
        newlog=1;
    end

    STEP=uint64(mem*(1024^3)/(ncols*8));
    m.logparams_all_trun(maxrow,ncols)=0;
    disp("starting truncation")
    while i+STEP<=maxrow
        try
            x=m.params_all(i:i+STEP,:);
            m.params_all_trun(i:i+STEP,:)=x;
            if newlog
                y=log10(x);
            else
                y=m.logparams_all(i:i+STEP,:);
            end
            m.logparams_all_trun(i:i+STEP,:)=y;
        catch
            disp("unable to load so trying to loop through each")
            for j=i:uint64(i+STEP)
                err=0;
                try
                    if newlog
                        x=m.params_all(j,:);
                        m.params_all_trun(j,:)=x;
                        m.logparams_all_trun(j,:)=log10(x);
                    else
                        try
                            x=m.params_all(j,:);
                        catch
                            err=1;
                        end
                        if err
                            y=m.logparams_all(j,:);
                            m.logparams_all_trun(j,:)=y;
                            m.params_all_trun(j,:)=10.^(y);
                        else
                            m.params_all_trun(j,:)=x;
                            m.logparams_all_trun(j,:)=log10(x);
                        end
                    end
                catch
                    m.params_all_trun(j,:)=NaN;
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
        x=m.params_all(ind,:);
        m.params_all_trun(ind,:)=x;
        if newlog
            y=log10(x);
        else
            y=m.logparams_all(ind,:);
        end
        m.logparams_all_trun(ind,:)=y;
    catch
        for j=ind
            err=0;
            try
                if newlog
                    x=m.params_all(j,:);
                    m.params_all_trun(j,:)=x;
                    m.logparams_all_trun(j,:)=log10(x);
                else
                    try
                        x=m.params_all(j,:);
                    catch
                        err=1;
                    end
                    if err
                        y=m.logparams_all(j,:);
                        m.logparams_all_trun(j,:)=y;
                        m.params_all_trun(j,:)=10.^(y);
                    else
                        m.params_all_trun(j,:)=x;
                        m.logparams_all_trun(j,:)=log10(x);
                    end
                end
            catch
                m.params_all_trun(j,:)=NaN;
                m.logparams_all_trun(j,:)=NaN;
                corruptindex=[corruptindex; j];
            end
        end
    end
    m.corruptindex=corruptindex;

    disp("params_all_trun and logparams_all_trun successfully made")

    if isempty(who(m,'params_out'))||length(m.params_out)==1
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
                lowerbound=uint64(minrow+i-1)
                %disp(lowerbound)
                upperbound=uint64(minrow+STEP+i-1)
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
        disp("taking log10 of params_out")
        m.parameters_all=log10(m.params_out);
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
    
    
    m.params_all=[];

    if ~isempty(who(m,'logparams_all'))
        m.logparams_all=[];
    end

    disp("removed params_all and logparams_all")

    m.updatedone=1;

end