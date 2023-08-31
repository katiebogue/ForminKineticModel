classdef Options <handle
    %OPTIONS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        python_path string
        scaling dictionary
        kpoly_type
    end

    methods
        function obj = Options(python_path)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.python_path=python_path;
        end

        function set_scaling(obj,step,scale_val,scale_type)
            arguments
                obj
                step (1,1) string {mustBeMember(step,{"kcap","kdel","rcap","rdel","krel"})}
                scale_val (1,1) string {mustBeMember(scale_val,{"dist_FH2","dist_NT", "size"})}
                scale_type (1,1) string {mustBeMember(scale_type,{"linear","exp"})}
            end
            if scale_type =="linear"
                x=@(PRM) PRM.(scale_val);
            elseif scale_type=="exp"
                x=@(PRM) exp(-1*PRM.(scale_val));
            end
            obj.scaling{(step)}=@(PRM) x(PRM);
        end
    end

end