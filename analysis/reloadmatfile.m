function reloadmatfile(matfileloc,mem,replaceTF)
    arguments
        matfileloc
        mem=15.5 % GB of RAM
        replaceTF=1
    end
    disp(matfileloc)
    disp(mem)
    disp(replaceTF)
    try
        if mem>6
            load(matfileloc)
            disp("successful load")
            save(matfileloc)
        else
            error("memory not high enough to attempt load")
        end
    catch
        disp("load failed, trying variable by variable")
        %save(strcat(matfileloc,'temp.mat'),'-v7.3')
        varlist=whos('-file', matfileloc);
        largein=[];
        for i=1:length(varlist)
            var=varlist(i).name;
            if var=="params_all_trun" || var=="params_out" || var=="params_all"
                % get rid of these parameters for older versions of the
                % code
            elseif varlist(i).bytes<=mem*(1024^3)
                load(matfileloc,var)
                %save(strcat(matfileloc,'temp.mat'),"var",'-append')
                fprintf("loaded var: %s\n",var)
                %clear("var")
            elseif var=="logparams_all_trun" || var=="parameters_all"
                largein=[largein, i];
            else
                disp(var)
                error("missing variable")
            end
        end
        save(strcat(matfileloc,'temp.mat'),'-v7.3')
        disp("saved workspace")
        for i=largein
            var=varlist(i).name;
            mnew=matfile(strcat(matfileloc,'temp.mat'),"Writable",true);
            mold=matfile(matfileloc,"Writable",false);
            sz=uint64(varlist(i).size);
            maxrows=mem*(1024^3)/(sz(2)*8);
            j=uint64(sz(1));
            while j-maxrows>0
                if var=="params_all_trun"
                    mnew.params_all_trun(j-maxrows:j,1:sz(2))=mold.params_all_trun(j-maxrows:j,1:sz(2));
                elseif var=="logparams_all_trun"
                    mnew.logparams_all_trun(j-maxrows:j,1:sz(2))=mold.logparams_all_trun(j-maxrows:j,1:sz(2));
                elseif var=="parameters_all"
                    mnew.parameters_all(j-maxrows:j,1:sz(2))=mold.parameters_all(j-maxrows:j,1:sz(2));
                end
                j=j-maxrows-1;
            end
            if var=="params_all_trun"
                mnew.params_all_trun(1:j,1:sz(2))=mold.params_all_trun(1:j,1:sz(2));
            elseif var=="logparams_all_trun"
                mnew.logparams_all_trun(1:j,1:sz(2))=mold.logparams_all_trun(1:j,1:sz(2));
            elseif var=="parameters_all"
                mnew.parameters_all(1:j,1:sz(2))=mold.parameters_all(1:j,1:sz(2));
            end
            clear mnew
            clear mold
            fprintf("saved looping var: %s\n",var)
        end
        if replaceTF
            movefile(strcat(matfileloc,'temp.mat'),matfileloc)
            delete(strcat(matfileloc,'temp.mat'))
        end
    end

end