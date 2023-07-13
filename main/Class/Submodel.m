classdef Submodel
    %SUBMODEL Contains values for single, double, and dimer submodels
    %   construct with: submodel(single,double_a,double_b,dimer_a,dimer_b)
    %   See also FILAMENT

    properties
        single double % contains value for single filament submodel
        double filament % value for double filament submodel (filament class)
        dimer  filament % value for dimer filament submodel (filament class)
    end
    properties (Dependent)
        ratio filament % ratio of dimer/ double values (filament class)
    end

    methods
        function obj = Submodel(single,double_a,dimer_a,double_b,dimer_b)
            %SUBMODEL Construct an instance of Submodel
            %   assign single to single
            %   assign double to filament of double_a, double_b (or both
            %   double_a if no b provided)
            %   assign dimer to filament of dimer_a, dimer_b (or both
            %   dimer_a if no b provided)
            if nargin==5
                obj.single=single;
                obj.double=filament(double_a,double_b);
                obj.dimer=filament(dimer_a,dimer_b);
            elseif nargin>0
                obj.single=single;
                obj.double=filament(double_a,double_a);
                obj.dimer=filament(dimer_a,dimer_b);
            end
        end

        function value = get.ratio(obj)
            value=obj.dimer/obj.double;
        end
    end
end