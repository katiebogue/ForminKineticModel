classdef Formin < handle
    %FORMIN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name string
        gating double=1
        length double
        c_PA double=1 % concentration of profilin-actin | μM
        params Params
        opts Options
    end

    properties (SetAccess=protected)
        sequence string="-1" % sequence that has been provided (not from UNIPROT)
        uniprotID string="-1"
        opts_over % overriden sequence options
        PRMList (1,:) PRM
        Ploc (1,:) double % locations of all Ps in sequence (relative to FH2)
    end

    properties(Dependent)
        kcap FilType
        kdel FilType
        kpoly FilType
        rcap FilType
        rdel FilType
        krel FilType
        extrapolation FilType
        PRMCount double
        meanPRMsize double
        numPs double % number of Ps IN PRMs (meanPRMsize x PRMCount)
        lookup Lookuptable
    end

    methods
        function obj = Formin(name,params,opts,NameValueArgs)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                name string
                params Params
                opts Options
                NameValueArgs.c_PA double
                NameValueArgs.gating double
                NameValueArgs.sequence string
                NameValueArgs.uniprotID string
            end
            if nargin>0
                obj.name=name;
                obj.params=params;
                obj.opts=opts;
                if isfield(NameValueArgs,"c_PA")
                    obj.c_PA=NameValueArgs.c_PA;
                end
                if isfield(NameValueArgs,"gating")
                    obj.gating=NameValueArgs.gating;
                end
                if isfield(NameValueArgs,"sequence") && isfield(NameValueArgs,"uniprotID")
                    error("cannot provide both an input sequence and a uniprotID")
                elseif isfield(NameValueArgs,"sequence")
                    obj.sequence=NameValueArgs.sequence;
                    obj.uniprotID="-1";
                elseif isfield(NameValueArgs,"uniprotID")
                    obj.uniprotID=NameValueArgs.uniprotID;
                    obj.sequence="-1";
                end
                addlistener(opts,'min_lenPRM','PostSet',@obj.update_FH1);
                addlistener(opts,'nInt','PostSet',@obj.update_FH1);
                addlistener(opts,'max_lenInt','PostSet',@obj.update_FH1);
                addlistener(opts,'min_nP','PostSet',@obj.update_FH1);
                addlistener(opts,'NTopt','PostSet',@obj.update_FH1);
                addlistener(opts,'CTopt','PostSet',@obj.update_FH1);
                addlistener(opts,'PRMscanopt','PostSet',@obj.update_FH1);
                addlistener(opts,'lookup','PostSet',@obj.update_FH1);
                obj.update_FH1();
            end
        end

        function value = get.extrapolation(obj)
            exsingle=isextrapolated("single");
            exdouble=isextrapolated("double");
            exdimer=isextrapolated("dimer");
            value=FilType(exsingle,exdouble,exdimer);

            function out = isextrapolated(type)
                if obj.length > obj.lookup.max_length.(type)
                    out=true;
                elseif ~isempty(obj.lookup.missingNs.(type))
                    out= ismember(strcat("N",num2str(obj.length)),obj.lookup.missingNs.(type));
                else
                    out=false;
                end
            end
        end


        function value=get.kcap(obj)
            value=get_rate(obj,"kcap");
        end

        function value=get.kdel(obj)
            value=get_rate(obj,"kdel");
        end

        function value=get.rcap(obj)
            value=get_rate(obj,"rcap");
        end

        function value=get.rdel(obj)
            value=get_rate(obj,"rdel");
        end

        function value=get.krel(obj)
            value=get_rate(obj,"krel");
        end

        function value=get.kpoly(obj)
            value=get_rate(obj,"kpoly");
        end

        function r=get_rate(obj,rate)
            rate_sum=FilType;
            for i=1:obj.PRMCount
                if i==1
                    rate_sum=obj.PRMList(i).(rate);
                else
                    rate_sum=rate_sum+obj.PRMList(i).(rate);
                end
            end
            if class(rate_sum)=="FilType"
                r=rate_sum.add_fils;
            else
                r=rate_sum;
            end
        end

        function value=get.PRMCount(obj)
            value=length(obj.PRMList);
        end

        function value=get.lookup(obj)
            value=obj.opts.lookup;
        end

        function value=get.meanPRMsize(obj)
            sum=0;
            for i=1:obj.PRMCount
                sum=sum+obj.PRMList(i).size;
            end
            value=sum/obj.PRMCount;
        end

        function value=get.numPs(obj)
            value=obj.PRMCount * obj.meanPRMsize;
        end
        
        function update_FH1(obj,src,evnt, NameValueArgs)
            arguments
                obj Formin
                src=1
                evnt=1
                NameValueArgs.sequence string
                NameValueArgs.uniprotID string
                NameValueArgs.min_lenPRM double % minimum PRM length (including interruptions)
                NameValueArgs.nInt double % number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions)
                NameValueArgs.max_lenInt double % maximum interruption length
                NameValueArgs.min_nP double % minimum number of Ps 
                NameValueArgs.NTopt double % FH1 NT definition options 
                NameValueArgs.CTopt double % FH1 CT definition options
                NameValueArgs.PRMscanopt double % PRM scan options
            end
            
            if isfield(NameValueArgs,"sequence") && isfield(NameValueArgs,"uniprotID")
                error("cannot provide both an input sequence and a uniprotID")
            elseif isfield(NameValueArgs,"sequence")
                obj.sequence=NameValueArgs.sequence;
                obj.uniprotID="-1";
            elseif isfield(NameValueArgs,"uniprotID")
                obj.uniprotID=NameValueArgs.uniprotID;
                obj.sequence="-1";
            elseif obj.sequence=="-1" && obj.uniprotID=="-1"
                return
            end
            
            min_lenPRM=getopt("min_lenPRM");
            nInt=getopt("nInt");
            max_lenInt=getopt("max_lenInt");
            min_nP=getopt("min_nP");
            NTopt=getopt("NTopt");
            CTopt=getopt("CTopt");
            PRMscanopt=getopt("PRMscanopt");

            output=runpyfunction(obj.opts,{'bioservices'},'get_formin_info','main',{obj.uniprotID, obj.sequence, min_lenPRM, nInt, max_lenInt, min_nP, NTopt, CTopt, PRMscanopt});
            obj.length=double(output(1));
            obj.Ploc=double(output{4});
            
            pp_index_vec=double(output{2});
            pp_length_vec=double(output{3});
            pp_index_vec_start=double(output{5});

            PRMs=PRM.empty;

            for i=1:length(output{2})
                dist_FH2=pp_index_vec(i);
                dist_NT=obj.length-pp_index_vec(i)+1;
                size=pp_length_vec(i);
                dist_FH2_start=pp_index_vec_start(i);
                PRMi=PRM(obj,dist_FH2,dist_NT,size,dist_FH2_start);
                PRMs(1,i)=PRMi;
            end

            obj.PRMList=PRMs;
            
            function output=getopt(option)
                if isfield(NameValueArgs,option)
                    obj.opts_over.(option)=NameValueArgs.(option);
                    output=NameValueArgs.(option);
                elseif isfield(obj.opts_over,option)
                    output=obj.opts_over.(option);
                else
                    output=obj.opts.(option);
                end
            end
        end

        function fig=makeSchematic(obj)
            fig=filamentSchematic(obj);
        end

        function fig=makeKpolyBar(obj)
            fig=forminkpolybar(obj);
        end
    end
end