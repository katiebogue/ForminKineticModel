% GENERATE_SWEEP_DATA2 creates and saves table with rates for formins of
% various configurations across random parameter value sweeps

new_opts=false;
    %% generate table of possible opts
    if new_opts
        opts_tab=add_opts(opts,1);
        opts_tab = table('Size',[3000 length(opts_tab.Properties.VariableNames)],'VariableNames',opts_tab.Properties.VariableNames,'VariableTypes',varfun(@class,opts_tab,'OutputFormat','cell'));
       for i=1:1000
           opts.cleareq;
        r=(20--4)*rand(1,5)+-4;
        opts.k_cap=10^(r(1));
        opts.k_del=10^(r(2));
        opts.r_cap=10^(r(3));
        opts.r_del=10^(r(4));
        opts.k_rel=10^(r(5));
        opts.kpoly_type="4st";
        opts.set_equation(1);
    
        opts_tab(i,:)=add_opts(opts,i);
    
        opts.set_equation(0,"krel",{"size","negexp"})
        opts_tab(1000+i,:)=add_opts(opts,1000+i);
    
        opts.kpoly_type="3st";
        opts_tab(2000+i,:)=add_opts(opts,2000+i);
       end
    end
   %% sweep through other possibilities
   for i=3004:3:3007
       lt=opts.lookup;
       clear opts
       opts=Options(lt,"","4st","",...   
        opts_tab{i,"k_cap"},...          % k_cap
        opts_tab{i,"k_del"},... % k_del
        opts_tab{i,"r_cap"},...    % r_cap
        opts_tab{i,"r_del"},...    % r_del
        opts_tab{i,"k_rel"});    % k_rel

       % opts=Options(lt,"","3st","",...   
       %  opts_tab{i,"k_cap"},...          % k_cap
       %  opts_tab{i,"k_del"},... % k_del
       %  opts_tab{i,"r_cap"},...    % r_cap
       %  0,...    % r_del
       %  0);    % k_rel
       opts.set_equation(1);
       clear lt
       % if i>2000
       %     opts.kpoly_type="3st";
       % else
       %  opts.set_equation(0,"krel",{"size","negexp"})
       % end
       opts.set_equation(0,"krel",{"size","negexp"})
       tab=zeros(1500,16);
       index=0;
       for FH1_length=10:20:600
        for PRM_size=[1,3,4,5,6,7,9,10,15,20]
            for PRM_loc=[ceil(PRM_size/2),ceil(FH1_length*0.25),ceil(FH1_length/2),ceil(FH1_length*0.75),FH1_length]
                formin1=Formin("formin1",opts,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
                kpoly=formin1.kpoly;
                kcap=formin1.kcap;
                kdel=formin1.kdel;
                rcap=formin1.rcap;
                rdel=formin1.rdel;
                krel=formin1.krel;
                % if i>2000
                %     rdel=0;
                %     krel=0;
                % else
                %     rdel=formin1.rdel;
                %     krel=formin1.krel;
                % end
                index=index+1;
                tab(index,:)=[i,FH1_length,PRM_loc,PRM_size,kpoly.double,kpoly.dimer,kpoly.ratio,kcap.double,kcap.dimer,kcap.ratio,kdel.double,kdel.dimer,kdel.ratio,rcap,rdel,krel];
                %tab(index,:)=[i,FH1_length,PRM_loc,PRM_size,kpoly.double,kpoly.dimer,kpoly.ratio,kcap.double,kcap.dimer,kcap.ratio,kdel.double,kdel.dimer,kdel.ratio,rcap,0,0];

            end
        end
       end
       disp(i)
       save(strcat("table_",num2str(i),".mat"),"tab",'-mat')
       clear tab
       clear formin1
       clear kcap
       clear kpoly
       clear kdel
   end

   %% 

function tab=add_opts(opts,index)
    tab=table;
    tab.id=index;
    tab.kpoly_type=opts.kpoly_type;
    
    params={"kcap_eq","kdel_eq","rcap_eq","rdel_eq","krel_eq"};
    for i=1:length(params)
        if isfield(opts.equationstext,params{i})
            tab.(params{i})=opts.equationstext.(params{i});
        else
            tab.(params{i})="n/a";
        end
    end
    
    params={"k_cap","k_del","r_cap","r_del","k_rel"};
    for i=1:length(params)
        tab.(params{i})=opts.(params{i});
    end

end

function savetable(tab,num)
    save(strcat("table_",num2str(num),".mat"),"tab",'-mat')
end
