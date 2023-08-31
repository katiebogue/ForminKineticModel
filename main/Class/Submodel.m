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
            elseif nargin==1
                obj.single=single;
                obj.double=single;
                obj.dimer=single;
            elseif nargin>0
                obj.single=single;
                obj.double=double_a;
                obj.dimer=dimer_a;
            end
        end

        function value = get.ratio(obj)
            value=obj.dimer/obj.double;
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
            singler=exp(obj.single);
            doubler=exp(obj.double);
            dimerr=exp(obj.dimer);
            r=Submodel(singler,doubler,dimerr);
        end


    end

    methods(Access = protected)
        function r = use_operator(obj1,obj2,operator)
            arguments
                obj1 {mustBeA(obj1,["double","Submodel"])}
                obj2 {mustBeA(obj2,["double","Submodel"])}
                operator function_handle
            end
            if class(obj1)=="double"
                singler=operator(obj1,obj2.single);
                doubler=operator(obj1,obj2.double);
                dimerr=operator(obj1,obj2.dimer);
            elseif class(obj2)=="Submodel"
                singler=operator(obj1.single,obj2.single);
                doubler=operator(obj1.double,obj2.double);
                dimerr=operator(obj1.dimer,obj2.dimer);
            elseif class(obj2)=="double"
                singler=operator(obj1.single,obj2);
                doubler=operator(obj1.double,obj2);
                dimerr=operator(obj1.dimer,obj2);
            end
            r=Submodel(singler,doubler,dimerr);
        end
    end

    methods(Access={?Formin})
        function obj=add_fils(obj)
            obj.double=obj.double.a+obj.double.b;
            obj.dimer=obj.dimer.a+obj.dimer.b;
        end
    end
end