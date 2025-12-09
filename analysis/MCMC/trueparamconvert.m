function params_all = trueparamconvert(mcmcfile)
%TRUEPARAMCONVERT convert nondimensionalized parameters into orignal scales
%
%   params_all = TRUEPARAMCONVERT(mcmcfile) convert nondimensionalized
%   parameters in the specified file into the orignal parameter values
%
%   Inputs:
%       mcmcfile : .mat file result from mcmcparamfit, must have the
%       following parameters: nparams, nkpolyparams, nsigma, type, nondim,
%       parameters_all, rates, divdatapoint, divkpoly 
%       
%   Outputs:
%       params_all : matrix with orignally scaled parameters
%
%   Loads values from .mat file, computes scaled kpoly values for each
%   parameter, and then converts the parameters into regular dimensions.
%
%   If the fit did not use nondimensionalization (nondim=0), then
%   params_all simply returns the parameters_all variable in the .mat file.
%
% See also MCMC_PRED.
    load(mcmcfile,'nparams')
    load(mcmcfile,'nkpolyparams')
    load(mcmcfile,'nsigma')
    load(mcmcfile,'type')
    load(mcmcfile,'nondim')
    load(mcmcfile,'parameters_all')
    load(mcmcfile,'rates')
    load(mcmcfile,'divdatapoint')
    load(mcmcfile,'divkpoly')
    load(mcmcfile,'prfit')
    load(mcmcfile,'prcalc')
    load(mcmcfile,'xloc')
    load(mcmcfile,'yloc')
    load(mcmcfile,'prmlocs')
    load(mcmcfile,'fh1lengths')

    if nondim==0
        params_all=parameters_all;
        return
    else
        params_all=zeros(length(parameters_all),nparams+1);
    end

    rates=struct2table(rates);
    rates=rates(divdatapoint,:);
    prmlocs=prmlocs{divdatapoint};
    fh1lengths=fh1lengths{divdatapoint};

    if prcalc
        if prfit
            
        else
            x1= xloc;
            y1= yloc;
            prdobs1=pr(prmlocs,fh1lengths,35.5,1,x1,y1,"double",1);
            prdims1=pr(prmlocs,fh1lengths,35.5,1,x1,y1,"dimer",1);
            for k=1:length(rates.k_delbase)
                prdob1=prdobs1;
                prdim1=prdims1;
                vals1=rates.k_delbase{k};
                for j=1:size(vals1,1)
                    vals1(j,1)=vals1(j,1)*prdob1(j);
                    vals1(j,2)=vals1(j,2)*prdob1(j);
                    vals1(j,3)=vals1(j,3)*prdob1(j);
                    vals1(j,4)=vals1(j,4)*prdim1(j);
                    vals1(j,5)=vals1(j,5)*prdim1(j);
                end
                rates.k_delbase{k}=vals1;
            end
        end
    end

    for i=1:length(parameters_all)
        kpoly_scale=calckpolyscale(parameters_all(i,:), prcalc, prfit, fh1lengths, prmlocs, rates);
        params_all(i,:)=gettrueparams(parameters_all(i,:),kpoly_scale,divkpoly,nkpolyparams, prfit);
        % if ~isreal(params_all)
        %     fh1lengths=fh1lengths{divdatapoint};
        % end
    end

    function kpoly_scale=calckpolyscale(params, prcalc, prfit, fh1lengths, prmlocs, rates)
        % calculate per PRM rates
        % params = alpha_del, beta_cap, rcapp_exp, (gamma_del, tau_rel)
 
        % calulcates kpolys for input parameters
        if prcalc
            if prfit
                x=params(end-1);
                y=params(end);
                prdobs=pr(prmlocs,fh1lengths,35.5,1,x,y,"double",1);
                prdims=pr(prmlocs,fh1lengths,35.5,1,x,y,"dimer",1);
                for k=1:length(rates.k_delbase)
                    prdob=prdobs;
                    prdim=prdims;
                    vals=rates.k_delbase{k};
                    for j=1:size(vals,1)
                        vals(j,1)=vals(j,1)*prdob(j);
                        vals(j,2)=vals(j,2)*prdob(j);
                        vals(j,3)=vals(j,3)*prdob(j);
                        vals(j,4)=vals(j,4)*prdim(j);
                        vals(j,5)=vals(j,5)*prdim(j);
                    end
                    rates.k_delbase{k}=vals;
                end
            end
        end
    
        % calculate per PRM rates
        
            % params = alpha_del, deta_cap, rcapp_exp, (gamma_del, tau_rel)
            % kcaps=cellfun(@(x) x, rates.k_capbase,'UniformOutput',false); 
            % kdels=cellfun(@(x) x.*10.^params(1), rates.k_delbase,'UniformOutput',false); 
            % rcaps=cellfun(@(x) ((x).^params(3)).*10.^params(2), rates.r_capbase,'UniformOutput',false);

            kcaps = cell(size(rates.k_capbase));
            for n = 1:numel(rates.k_capbase)
                kcaps{n} = rates.k_capbase{n};
            end
    
            kdels = cell(size(rates.k_delbase));
            for n = 1:numel(rates.k_delbase)
                kdels{n} = rates.k_delbase{n} .* 10.^params(1);
            end
    
            rcaps= cell(size(rates.r_capbase));
            for n = 1:numel(rates.r_capbase)
                rcaps{n} = ((rates.r_capbase{n}).^params(3)).*10.^params(2);
            end
        if type=="4st"
                %rdels=cellfun(@(x) x.*10.^params(4), rates.r_delbase,'UniformOutput',false);
                rdels = cell(size(rates.r_delbase));
                for n = 1:numel(rates.r_delbase)
                    rdels{n} = rates.r_delbase{n} .*10.^params(4);
                end

                %krels=cellfun(@(x) x.*10.^params(5), rates.k_relbase,'UniformOutput',false);
                krels = cell(size(rates.k_relbase));
                for n = 1:numel(rates.k_relbase)
                    kdels{n} = rates.k_relbase{n} .*10.^params(5);
                end
            
            %kpolys=cellfun(@(kcap,kdel,rcap,rdel,krel) 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel))),kcaps,kdels,rcaps,rdels,krels,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
            kpolys = cell(size(kcaps));
            for n = 1:numel(kcaps)
                kcap = kcaps{n};
                kdel = kdels{n};
                rcap = rcaps{n};
                rdel = rdels{n};
                krel = krels{n};
                kpolys{n} = 1./((1./krel) + ((rdel + krel)./(kdel .* krel)) + (((rcap .* rdel) + (rcap .* krel) + (kdel .* krel))./(kcap .* kdel .* krel)));
            end
        elseif type=="3st"
            %rdels=kcaps;
            %krels=kcaps;
            %kpolys=cellfun(@(kcap,kdel,rcap) 1./((1./kdel) + ((kdel + rcap)./(kdel.*kcap))),kcaps,kdels,rcaps,'UniformOutput',false); % using formin inputs, calculate double and dimer for all formins
            kpolys = cell(size(kcaps));
            for n = 1:numel(kcaps)
                kcap = kcaps{n};
                kdel = kdels{n};
                rcap = rcaps{n};
                kpolys{n} = 1 ./ ((1 ./ kdel) + ((kdel + rcap) ./ (kdel .* kcap)));
            end
        end

           
        
    
        for j=1:length(kpolys)
            PRMsum=sum(kpolys{j},1); % sum up PRMs
            kpolys{j}=[PRMsum(1),sum(PRMsum(2:3)),sum(PRMsum(4:5))]; % Sum up filaments
        end

        kpoly_scale=kpolys{1}(2);
        % if kpoly_scale<0
        %     return
        % end
        if kpoly_scale<9^-100
            return
        end
    end
end

function trueparams=gettrueparams(params,alphakp,kp,nparams, prfit)
    %must be 3 state method
    kcap=kp/alphakp;
    trueparams=log10((10.^params)*kcap);
    trueparams(3)=params(3); %rcap_exp

    if prfit
        trueparams(nparams-1:nparams)=params(nparams-1:nparams); %delivery locations
    end

    i=nparams;
    while i<length(trueparams)
        % dont change sigma values
        i=i+1;
        trueparams(i)=params(i);
    end

    trueparams=[log10(kcap), trueparams];
    % if ~isreal(trueparams)
    %     return 
    % end
    if isinf(trueparams(2))
        return
    end
end