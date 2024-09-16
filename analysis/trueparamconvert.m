function params_all = trueparamconvert(mcmcfile)
    load(mcmcfile,'nparams')
    load(mcmcfile,'nkpolyparams')
    load(mcmcfile,'nsigma')
    load(mcmcfile,'type')
    load(mcmcfile,'nondim')
    load(mcmcfile,'parameters_all')
    load(mcmcfile,'rates')
    load(mcmcfile,'divdatapoint')
    load(mcmcfile,'divkpoly')

    if nondim==0
        params_all=parameters_all;
        return
    else
        params_all=zeros(length(parameters_all),nparams+1);
    end

    rates=struct2table(rates);
    rates=rates(divdatapoint,:);

    for i=1:length(parameters_all)
        kpoly_scale=calckpolyscale(parameters_all(i,:));
        params_all(i,:)=gettrueparams(parameters_all(i,:),kpoly_scale,divkpoly,nkpolyparams);
    end

    function kpoly_scale=calckpolyscale(params)
        % calculate per PRM rates
        % params = alpha_del, deta_cap, rcapp_exp, (gamma_del, tau_rel)
        kcaps=cellfun(@(x) x, rates.k_capbase,'UniformOutput',false); 
        kdels=cellfun(@(x) x.*params(1), rates.k_delbase,'UniformOutput',false); 
        rcaps=cellfun(@(x) ((x).^params(3)).*params(2), rates.r_capbase,'UniformOutput',false); 
        if type=="4st"
            rdels=cellfun(@(x) x.*params(4), rates.r_delbase,'UniformOutput',false);
            krels=cellfun(@(x) x.*params(5), rates.k_relbase,'UniformOutput',false);
            kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
        elseif type=="3st"
            rdels=kcaps;
            krels=kcaps;
            kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
        end
        
    
        for j=1:length(kpolys)
            PRMsum=sum(kpolys{j},1); % sum up PRMs
            kpolys{j}=[PRMsum(1),sum(PRMsum(2:3)),sum(PRMsum(4:5))]; % Sum up filaments
        end

        kpoly_scale=kpolys{1}(2);
    end

end


function trueparams=gettrueparams(params,alphakp,kp,nparams)
    %must be 3 state method
    kcap=kp/alphakp;
    trueparams=log10((10.^params)*kcap);
    trueparams(3)=params(3); %rcap_exp

    i=nparams;
    while i<length(trueparams)
        % dont change sigma values
        i=i+1;
        trueparams(i)=params(i);
    end

    trueparams=[log10(kcap), trueparams];
end