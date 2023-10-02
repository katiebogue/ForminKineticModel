classdef Experiment 
    %EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ForminList (:,:) Formin
        opts Options
        params Params
    end

    methods
        function obj = Experiment(options,params,file,type,cPA)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options Options
                params Params
                file string % .txt file, comma separated
                type string {mustBeMember(type,{'seq','uniprot'})}
                cPA double =2.5 % [profilin-actin] to assign to all formins
            end
            if nargin>0
                obj.opts=options;
                obj.params=params;
                forminlist = char(importdata(file)); 
                forminlist = strsplit(forminlist);
                obj.ForminList=Formin.empty((length(forminlist)/2),0);
                for i = 1:length(forminlist)/2
                    forminname = convertCharsToStrings(forminlist(2*i -1));   %takes the names (every other string)
                    input = convertCharsToStrings(forminlist(2*i));  
                    if type=="seq"
                        tempFormin=Formin(forminname,obj.params,obj.opts,c_PA=cPA,sequence=input);
                    elseif type=="uniprot"
                        tempFormin=Formin(forminname,obj.params,obj.opts,c_PA=cPA,uniprotID=input);
                    end
                    obj.ForminList(i)=tempFormin;
                end
            end
        end

        function obj=add_formin(obj,formin)
            arguments
                obj Experiment
                formin Formin
            end
            obj.ForminList(1,end+1)=formin;
        end

        function obj= set.opts(obj,input)
            obj.opts=input;
            obj.syncforminsettings;
        end

        function obj= set.params(obj,input)
            obj.params=input;
            obj.syncforminsettings;
        end

        function syncforminsettings(obj)
            if isfield(obj,"ForminList")
                for i=1:length(obj.ForminList)
                    if obj.ForminList(i).opts~=obj.opts
                        obj.ForminList(i).opts=obj.opts;
                    end
                    if obj.ForminList(i).params~=obj.params
                        obj.ForminList(i).params=obj.params;
                    end
                end
            end
        end
    end
end