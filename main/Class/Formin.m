classdef Formin < handle
    %FORMIN Contains information and functions for a formin FH1 domain.
    %   
    %   Construction:
    %       obj = FORMIN(name, opts, NameValueArgs) 
    %
    %       NameValueArgs:
    %           c_PA     : (double) concentration of profilin actin
    %                     (NameValueArgs)
    %           gating   : (double) gating factor (NameValueArgs)
    %           sequence : (string) input sequence, can't have a uniprotID (NameValueArgs)
    %           uniprotID: (string) UNIPROT ID number, can't have an input sequence (NameValueArgs)
    %           length   : (double) FH1 length, number of amino acids (NameValueArgs)
    %           PRMloc   : (double) location of PRM, must be provided
    %                    if giving length and PRMsize (NameValueArgs)
    %           PRMsize  : (double) number of amino acids in PRM, must be provided
    %                    if giving length and PRMloc (NameValueArgs)
    %       
    %       3 sequence options:
    %           obj = FORMIN(name, opts, 'sequence',seq) construct a formin 
    %                 with the specified input sequence seq 
    %           
    %           obj = FORMIN(name, opts, 'uniprotID',ID) construct a formin
    %                 using the sequence of the Uniprot entry with the 
    %                 UniprotID ID
    %           
    %           obj = FORMIN(name, opts, 'length',L,'PRMloc',l,'PRMsize',n) 
    %                 construct a formin of length L with a single PRM n
    %                 amino acids long and l amino acids away from the FH2
    %                 domain
    %
    % See also PRM, EXPERIMENT, OPTIONS.

    properties
        name string % name/identifier for the formin
        gating double=1 % gating factor
        length double % number of amino acids in the FH1
        c_PA double=1 % concentration of profilin-actin | μM
        opts Options % Options object to use
    end

    properties (SetAccess=protected)
        sequence string="-1" % sequence that has been provided (not from UNIPROT)
        uniprotID string="-1" % UNIPROT id (if not using an input sequence)
        opts_over % list of overriden sequence options (used instead of opts)
        PRMList (1,:) PRM % PRM's in the FH1, children
        Ploc (1,:) double % locations of all Ps in sequence (relative to FH2)
        lastkpoly FilType % most recently computed kpoly values, can be accessed without recalculating things, but might not be the most updated
        NName string % FH1 length in string format with leading "N"
    end

    properties(Dependent)
        kcap FilType % (sum of PRMs) rate of PRM + profilin-actin binding (capture)| s^(-1)
        kdel FilType % (sum of PRMs) rate of barbed end + PRM-profilin-actin binding (delivery) | s^(-1)
        kpoly FilType % (sum of PRMs) overall rate of polymerization | s^(-1)
        rcap FilType % (sum of PRMs) rate of PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        rdel FilType % (sum of PRMs) rate of barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        krel FilType % (sum of PRMs) rate of PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
        extrapolation FilType % whether or not (true/false) extrapolated polymer stats were used to calculate kpolys for each FilType
        PRMCount double % number of PRMs
        meanPRMsize double % average size (total number of amino acids) of PRMs
        numPs double % number of Ps IN PRMs (meanPRMsize x PRMCount)
        lookup Lookuptable % corresponding lookuptable, from opts property
    end

    methods
        function obj = Formin(name,opts,NameValueArgs)
            %FORMIN Construct an instance of Formin
            %    
            %   obj = FORMIN(name, opts, NameValueArgs) construct a formin 
            %         with name and specified options
            %   Inputs: 
            %       name     : (string) name of formin
            %       opts     : options object to use
            %       c_PA     : (double) concentration of profilin actin
            %                   (NameValueArgs)
            %       gating   : (double) gating factor (NameValueArgs)
            %       sequence : (string) input sequence, can't have a uniprotID (NameValueArgs)
            %       uniprotID: (string) UNIPROT ID number, can't have an input sequence (NameValueArgs)
            %       length   : (double) FH1 length, number of amino acids (NameValueArgs)
            %       PRMloc   : (double) location of PRM, must be provided
            %                  if giving length and PRMsize (NameValueArgs)
            %       PRMsize  : (double) number of amino acids in PRM, must be provided
            %                  if giving length and PRMloc (NameValueArgs)
            %   
            %   Assigns input values to properties, adds listeners for options, and then runs update_FH1
            %   
            %   See also FORMIN, OPTIONS, FORMIN/UPDATE_FH1.
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
            % return true/false values for if extrapolation was used to
            % compute kpoly for each filament type
            exsingle=isextrapolated("single");
            exdouble=isextrapolated("double");
            exdimer=isextrapolated("dimer");
            value=FilType(exsingle,exdouble,exdimer);

            function out = isextrapolated(type)
                % determines if extrapolation was used for the input type
                % by checking if the formin's length is included in the
                % lookuptable
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
            % compute kcap rate from PRMs
            value=get_rate(obj,"kcap");
        end

        function value=get.kdel(obj)
            % compute kdel rate from PRMs
            value=get_rate(obj,"kdel");
        end

        function value=get.rcap(obj)
            % compute rcap rate from PRMs
            value=get_rate(obj,"rcap");
        end

        function value=get.rdel(obj)
            % compute rdel rate from PRMs
            value=get_rate(obj,"rdel");
        end

        function value=get.krel(obj)
            % compute krel rate from PRMs
            value=get_rate(obj,"krel");
        end

        function value=get.kpoly(obj)
            % compute overall kpoly rate from PRMs
            value=get_rate(obj,"kpoly");
            obj.lastkpoly=value; % set lastkpoly to the computed value
        end

        function r=get_rate(obj,rate)
            % GET_RATE compute rate from sum of rates of each PRM
            %
            %   r=FORMIN.GET_RATE(rate) computes rate of the formin
            %
            %   Adds up the rate values for each PRM.
            %   Since the rates for each PRM are dependent, the individual
            %   PRM rates are re-calculated every time this function runs.
            %
            % See also PRM.
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
            % determine the number of PRMs in PRMList
            value=length(obj.PRMList);
        end

        function value=get.lookup(obj)
            % retrieve the lookuptable object from opts
            value=obj.opts.lookup;
        end

        function value=get.meanPRMsize(obj)
            % calculate the mean PRM size of the PRMs in PRMList
            sum=0;
            for i=1:obj.PRMCount
                sum=sum+obj.PRMList(i).size;
            end
            value=sum/obj.PRMCount;
        end

        function value=get.numPs(obj)
            % calculate the total number of prolines in the FH1
            value=obj.PRMCount * obj.meanPRMsize;
        end
        
        function add_length(obj,added_length)
            %ADD_LENGTH increase the FH1 length but keep PRMs in place
            %
            %   FORMIN.ADD_LENGTH(x) increase the formin length by x
            %
            %   Inputs:
            %       added_length: (double) number of amino acids to
            %       increase the FH1 length by (can be negative)
            %
            %   See also FORMIN, PRM.
            obj.length=obj.length+added_length;
            obj.NName=strcat("N",num2str(obj.length));
            for i=1:length(obj.PRMList)
                PRM=obj.PRMList(1,i);
                PRM.dist_NT=PRM.dist_NT+added_length;
                %PRM.updateStats;
            end
        end
        function update_FH1(obj,src,evnt, NameValueArgs)
            %UPDATE_FH1 read in formin sequence and update PRMs accordingly
            %
            %   FORMIN.UPDATE_FH1 re-read in formin sequence according to
            %   properties of opts and sequence properties of the formin
            %
            %   Inputs:
            %       sequence   : (string) input sequence, can't have a uniprotID (NameValueArgs)
            %       uniprotID  : (string) UNIPROT ID number, can't have an input sequence (NameValueArgs)
            %       min_lenPRM : (double) minimum PRM length (including interruptions), use instead of opts property (NameValueArgs)
            %       nInt       : (double) number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions), use instead of opts property (NameValueArgs)
            %       max_lenInt : (double) maximum interruption length, use instead of opts property (NameValueArgs)
            %       min_nP     : (double) minimum number of Ps in PRM, use instead of opts property (NameValueArgs)
            %       NTopt      : (double) FH1 NT definition options , use instead of opts property (NameValueArgs)
            %                   1 - first instance of PRM with at least 4 Ps with max 1 interruption of length 1 (in sequence of at least 3 PRMs no father than 100 amino acids apart)
            %                   2 - first instance of PRM (as defined by args 3 and 4) (uin sequence of at least 3 PRMs no father than 100 amino acids apart)
            %                   3 - Uniprot defined FH1 start (or option 1 if no FH1)
            %                   4 - Uniprot defined FH1 start (or option 2 if no FH1)
            %                   5 - Start of sequence (for input sequences)
            %       CTopt      : (double) FH1 CT definition options, use instead of opts property (NameValueArgs)
            %                   1 - Uniprot defined FH2 start
            %                   2 - Uniprot defined FH1 end (or option 1 if no FH1)
            %                   3 - End of sequence (for input sequences)
            %       PRMscanopt : (double) PRM scan options, use instead of opts property (NameValueArgs)
            %                   1 - Search for PRMs starting from the FH2 (CT) (default)
            %                   2 - Search for PRMs starting from the FH1 NT
            %
            %       Uses runpyfunction to run get_formin_info.py
            %
            %   See also FORMIN, RUNPYFUNCTION.
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
                % retrieve application options value, from opts or from
                % input if overriden via NameValueArgs
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
            % MAKESCHEMATIC run filamentSchematic on the Formin
            %
            %   fig = FORMIN.MAKESCHEMATIC creates schematic of formin. 
            %
            %   Schematic is 2 lines representing the entire length of an fh1
            %   filament. 1 line has a colored dot for the center of each binding 
            %   site; 1 line has blue dots for each proline, with binding sites 
            %   highlighted in yellow.
            % 
            %   See also FILAMENTSCHEMATIC.
            fig=filamentSchematic(obj);
        end

        function fig=makeKpolyBar(obj,scale)
            % MAKEKPOLYBAR run forminkpolybar on the Formin; creates bargraphs of the rates of each step in and of 
            % polymerization for the formin.
            %   
            %   fig = FORMIN.MAKEKPOLYBAR creates figure with subplot for the 
            %   rate of each step calculated for formin. 
            %
            %   fig = FORMIN.MAKEKPOLYBAR('scale') creates figure with subplot for the
            %   rate of each step (scaled as specified) calculated for formin. 
            %
            %   Determines the steps from formin.opts.equations
            %   Bargraph is stacked, colored by PRM
            %
            %   Inputs:
            %       scale  : (string) kpoly axis scale (can be none, log2, log10, ln; default is none)
            %
            % See also FORMINKPOLYBAR, OPTIONS
            arguments
                obj Formin
                scale string="none"
            end
            fig=forminkpolybar(obj,scale);
        end

        function fig=formingraphic(obj,save)
            % FORMINGRAPHIC create figure with formin schematic and kpoly
            % bar graphs
            %
            %   fig = FORMIN.FORMINGRAPHIC creates formin schematic and kpoly
            %   bar graphs. 
            %   
            %   fig = FORMIN.FORMINGRAPHIC(true) creates and saves formin schematic and kpoly
            %   bar graphs.
            %
            %   Runs formin.makeSchematic and formin.makeKpolyBar and adds
            %   both outputs to a tiled layout
            %
            %   Inputs:
            %       save : (logical) whether or not to save the figure (default is false)
            %
            %   Calls figuresave
            %
            % See also FORMIN/MAKESCHEMATIC, FORMIN/MAKEKPOLYBAR,
            % FIGURESAVE.
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
            %LOADOBJ run tasks on load of a Formin object
            %
            % On load, add new listeners for the formin's options object to trigger update_formin.
            %
            % See also FORMIN/UPDATE_FH1, OPTIONS.
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