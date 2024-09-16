classdef Formin < handle
    %FORMIN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name string
        gating double=1
        length double
        c_PA double=1 % concentration of profilin-actin | μM
        opts Options
    end

    properties (SetAccess=protected)
        sequence string="-1" % sequence that has been provided (not from UNIPROT)
        uniprotID string="-1"
        opts_over % overriden sequence options
        PRMList (1,:) PRM
        Ploc (1,:) double % locations of all Ps in sequence (relative to FH2)
        lastkpoly FilType
        NName string
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
        function obj = Formin(name,opts,NameValueArgs)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                name string
                opts Options
                NameValueArgs.c_PA double
                NameValueArgs.gating double
                NameValueArgs.sequence string
                NameValueArgs.uniprotID string
                NameValueArgs.length double
                NameValueArgs.PRMloc double
                NameValueArgs.PRMsize double
               
            end
            if nargin>0
                obj.name=name;
                obj.opts=opts;
                if isfield(NameValueArgs,"c_PA")
                    obj.c_PA=NameValueArgs.c_PA;
                end
                if isfield(NameValueArgs,"gating")
                    obj.gating=NameValueArgs.gating;
                end

                if isfield(NameValueArgs,"length")
                    if isfield(NameValueArgs,"sequence") || isfield(NameValueArgs,"uniprotID")
                        error("must either provide input values or input sequence/uniprot, not both")
                    end
                    if ~isfield(NameValueArgs,"PRMloc")
                        error("must provide PRMloc if using input values")
                    elseif ~isfield(NameValueArgs,"PRMsize")
                         error("must provide PRMsize if using input values")
                    else
                        obj.length=NameValueArgs.length;
                        obj.NName=strcat("N",num2str(obj.length));
                        obj.PRMList=PRM.empty;
                        PRMstartloc=NameValueArgs.PRMloc-ceil(NameValueArgs.PRMsize/2)+1;
                        obj.PRMList(1,1)=PRM(obj,NameValueArgs.PRMloc,(obj.length-NameValueArgs.PRMloc),NameValueArgs.PRMsize,PRMstartloc);
                        obj.Ploc=PRMstartloc:(PRMstartloc+NameValueArgs.PRMsize);
                    end
                elseif isfield(NameValueArgs,"sequence") && isfield(NameValueArgs,"uniprotID")
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
            obj.lastkpoly=value;
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
        
        function add_length(obj,added_length)
            obj.length=obj.length+added_length;
            obj.NName=strcat("N",num2str(obj.length));
            for i=1:length(obj.PRMList)
                PRM=obj.PRMList(1,i);
                PRM.dist_NT=PRM.dist_NT+added_length;
                PRM.dist_FH2_start=PRM.dist_FH2_start+added_length;
                %PRM.updateStats;
            end
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
            obj.NName=strcat("N",num2str(obj.length));
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

        function fig=makeKpolyBar(obj,scale)
            arguments
                obj Formin
                scale string="none"
            end
            fig=forminkpolybar(obj,scale);
        end

        function fig=formingraphic(obj,save)
            arguments
                obj Formin
                save logical=false
            end
            f1=obj.makeSchematic;
            f2=obj.makeKpolyBar;
            fig=figure;
            tiles=tiledlayout(1,length(f2.Children)+1);
            ax1=f1.Children;
            ax1.Parent=tiles;
            ax1.Layout.Tile=1;
            j=1;
            for i=length(f2.Children):-1:1
                j=j+1;
                ax=f2.Children(i);
                ax.Parent=tiles;
                ax.Layout.Tile=j;
            end
            title(tiles,obj.name)
            close(f1)
            close(f2)

            if save
                figuresave(gcf,obj.opts,append(obj.name,'.fig'));
            end

        end
    end

    methods (Static)
        function obj = loadobj(s)
            obj=s;
            addlistener(obj.opts,'min_lenPRM','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'nInt','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'max_lenInt','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'min_nP','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'NTopt','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'CTopt','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'PRMscanopt','PostSet',@obj.update_FH1);
            addlistener(obj.opts,'lookup','PostSet',@obj.update_FH1);
        end
    end
end