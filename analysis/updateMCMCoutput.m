function updateMCMCoutput(matfileloc)

    m = matfile(matfileloc,'Writable',true);

    if ~isempty(who(m,'updatedone'))
        if m.updatedone
            return
        end
    end

    if ~isempty(who(m,'paccept_matrix'))
        m.paccept_matrix( ~any(m.paccept_matrix,2),:) = [];
    end
    if ~isempty(who(m,'ksvals'))
        m.ksvals( ~any(m.ksvals,2),:) = [];
    end
    if ~isempty(who(m,'paramHistCounts_matrices'))
        m.paramHistCounts_matrices_nts( ~any(m.paramHistCounts_matrices_nts,2)) = [];
        m.paramHistCounts_matrices(:,:,length(m.paramHistCounts_matrices_nts)+1:end)=[];
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

    m.params_all_trun(maxrow,ncols)=0;
    corruptindex=[];
    i=1;

    newlog=0;
    if isempty(who(m,'logparams_all'))
        newlog=1;
    end

    STEP=100000;
    m.logparams_all_trun(maxrow,ncols)=0;
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
            for j=i:(i+STEP)
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
    end

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
    % for j=maxrow:-1:i
    %         err=0;
    %         try
    %             if newlog
    %                 try
    %                     x=m.params_all(j,:);
    %                 catch
    %                     err=1;
    %                 end
    %                 if err
    %                     y=m.logparams_all(j,:);
    %                     m.logparams_all_trun(j,:)=y;
    %                     m.params_all_trun(j,:)=10.^(y);
    %                 else
    %                     m.params_all_trun(j,:)=x;
    %                     m.logparams_all_trun(j,:)=log10(x);
    %                 end
    % 
    %             else
    %                 x=m.params_all(j,:);
    %                 m.params_all_trun(j,:)=x;
    %                 m.logparams_all_trun(j,:)=log10(x);
    %             end
    %         catch
    %             m.params_all_trun(j,:)=NaN;
    %             m.logparams_all_trun(j,:)=NaN;
    %             corruptindex=[corruptindex; j];
    %         end
    % end
    m.corruptindex=corruptindex;

    

    if isempty(who(m,'params_out'))||length(m.params_out)==1
        m.parameters_all=m.logparams_all_trun(ceil(maxrow/3):maxrow,:);
    else
        m.parameters_all=log10(m.params_out);
    end
    
    nparams=m.nparams;
    m.containsInf=zeros(nparams,1);
    len=size(m,'parameters_all',1);
    for i=1:nparams
        if max(m.parameters_all(1:len,i))==Inf || min(m.parameters_all(1:len,i))==-Inf
            m.containsInf(i,1)=1;
        end
    end

    m.params_all=[];

    if ~isempty(who(m,'logparams_all'))
        m.logparams_all=[];
    end

    m.updatedone=1;

end