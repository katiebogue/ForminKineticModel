classdef Filament
    %FILAMENT Contains values for a and b filament
    % construct with: filament(a_val,b_val)

    properties
        a  % value for filament a
        b  % value for filament b
    end
    properties(Dependent)
        avg % average of a and b values
    end

    methods
        function obj = Filament(a_val,b_val)
            %FILAMENT Construct an instance of Filament
            %   assigns a_val to a and b_val to b
            if nargin>0
                obj.a = a_val;
                obj.b= b_val;
            end
        end

        function value = get.avg(obj)
            arr=[obj.a,obj.b];
            value=mean(arr);
        end

        function r = mrdivide(obj1,obj2)
            operator=@mrdivide;
            r=use_operator(obj1,obj2,operator);
        end

        function r = plus(obj1,obj2)
            operator=@plus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = minus(obj1,obj2)
            operator=@minus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = mtimes(obj1,obj2)
            operator=@mtimes;
            r=use_operator(obj1,obj2,operator);
        end

        function r = mpower(obj1,obj2)
            operator=@mpower;
            r=use_operator(obj1,obj2,operator);
        end

        function r = exp(obj)
            ar=exp(obj.a);
            br=exp(obj.b);
            r=Filament(ar,br);
        end
    end

    methods(Access = protected)
        function r = use_operator(obj1,obj2,operator)
            arguments
                obj1 {mustBeA(obj1,["double","Filament"])}
                obj2 {mustBeA(obj2,["double","Filament"])}
                operator function_handle
            end
            if class(obj1)=="double"
                ar=operator(obj1,obj2.a);
                br=operator(obj1,obj2.b);
            elseif class(obj2)=="Filament"
                ar=operator(obj1.a,obj2.a);
                br=operator(obj1.b,obj2.b);
            elseif class(obj2)=="double"
                ar=operator(obj1.a,obj2);
                br=operator(obj1.b,obj2);
            end
            r=Filament(ar,br);
        end
    end
end