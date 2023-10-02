classdef FilType
    %FILTYPE Contains values for single, double, and dimer filament types
    %   construct with: FilType(single,double_a,double_b,dimer_a,dimer_b)
    %   See also FILAMENT

    properties
        single % contains value for single filament type
        double % value for double filament type (filament class if 2)
        dimer  % value for dimer filament type (filament class if 2)
        intersectratio logical=false % if true, ratio will instead be overlapping value between double and dimer
    end
    properties (Dependent)
        ratio Filament % ratio of dimer/ double values (filament class)
    end

    methods
        function obj = FilType(single,double_a,dimer_a,double_b,dimer_b)
            %FILTYPE Construct an instance of FilType
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
            value=getratio(obj);
            function out= getratio(input)
                if class(input.dimer)=="struct"
                    fieldsdimer=fieldnames(input.dimer);
                    fieldsdouble=fieldnames(input.double);
                    fields=intersect(fieldsdimer,fieldsdouble);
                    out=struct;
                    for i=1:length(fields)
                        temp1=FilType(0,input.double.(fields{i}),input.dimer.(fields{i}));
                        temp=getratio(temp1);                                                                                                                                                                                                                                
                        out.(fields{i})=temp;
                    end
                elseif class(input.dimer)=="Filament"
                    outa=getratio(FilType(0,input.double.a,input.dimer.a));
                    outb=getratio(FilType(0,input.double.b,input.dimer.b));
                    out=Filament(outa,outb);
                elseif class(input.dimer)=="string"
                    out="n/a";
                else
                    if input.double==0
                        out=NaN;
                    elseif input.dimer==0
                        out=0;
                    else
                        if obj.intersectratio
                            out=intersect(input.dimer,input.double);
                        else
                            out=input.dimer./input.double;
                        end
                    end
                end
            end
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
            r=FilType(singler,doubler,dimerr);
        end


    end

    methods(Access = protected)
        function r = use_operator(obj1,obj2,operator)
            arguments
                obj1 {mustBeA(obj1,["double","FilType"])}
                obj2 {mustBeA(obj2,["double","FilType"])}
                operator function_handle
            end
            if class(obj1)=="double"
                singler=operator(obj1,obj2.single);
                doubler=operator(obj1,obj2.double);
                dimerr=operator(obj1,obj2.dimer);
            elseif class(obj2)=="FilType"
                singler=operator(obj1.single,obj2.single);
                doubler=operator(obj1.double,obj2.double);
                dimerr=operator(obj1.dimer,obj2.dimer);
            elseif class(obj2)=="double"
                singler=operator(obj1.single,obj2);
                doubler=operator(obj1.double,obj2);
                dimerr=operator(obj1.dimer,obj2);
            end
            r=FilType(singler,doubler,dimerr);
        end
    end

    methods(Access={?Formin})
        function obj=add_fils(obj)
            obj.double=obj.double.a+obj.double.b;
            obj.dimer=obj.dimer.a+obj.dimer.b;
        end
    end
end