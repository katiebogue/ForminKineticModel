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
            obj.a = a_val;
            obj.b= b_val;
        end

        function value = get.avg(obj)
            arr=[obj.a,obj.b];
            value=mean(arr);
        end

        function r = mrdivide(obj1,obj2)
            ar=obj1.a/obj2.a;
            br=obj1.b/obj2.b;
            r=filament(ar,br);
        end
    end
end