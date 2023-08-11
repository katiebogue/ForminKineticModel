classdef Submodel
    %SUBMODEL Contains values for single, double, and dimer submodels
    %   construct with: submodel(single,double_a,double_b,dimer_a,dimer_b)
    %   See also FILAMENT

    properties
        single % contains value for single filament submodel
        double % value for double filament submodel (filament class if 2)
        dimer  % value for dimer filament submodel (filament class if 2)
    end
    properties (Dependent)
        ratio Filament % ratio of dimer/ double values (filament class)
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
                obj.double=Filament(double_a,double_b);
                obj.dimer=Filament(dimer_a,dimer_b);
            elseif nargin>0
                obj.single=single;
                obj.double=double_a;
                obj.dimer=dimer_a;
            end
        end

        function value = get.ratio(obj)
            value=obj.dimer/obj.double;
        end
    end
end