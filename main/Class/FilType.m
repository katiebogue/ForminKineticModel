classdef FilType
%FILTYPE Contains values for single, double, and dimer filament types
    %
    %   Construction:
    %       val=FILTYPE(single,double_a,dimer_a,double_b,dimer_b)
    %
    %       val=FILTYPE(single,double_a,dimer_a)
    %
    %   Multiple operators have been overloaded for this class (see methods).
    %
    %   See also FILAMENT.

    properties
        single % contains value for single filament type
        double % value for double filament type (filament class if 2)
        dimer  % value for dimer filament type (filament class if 2)
        intersectratio logical=false % if true, ratio will instead be overlapping values between double and dimer
    end
    properties (Dependent)
        ratio Filament % ratio of dimer/ double values (filament class; Dependent)
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
                elseif length(input.double)>1 && length(input.double)==length(input.dimer)
                    if obj.intersectratio
                            out=intersect(input.dimer,input.double);
                    else
                        out=input.double;
                        for i=1:length(input.double)
                            fil=FilType(0,input.double(i),input.dimer(i));
                            out(i)=fil.ratio;
                        end
                    end
                else
                    if all(input.double==0) && all(input.dimer==0)
                        out=1;
                    elseif input.double==0
                        out=0;
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
                obj1 
                obj2 
                operator function_handle
            end
            class1=class(obj1);
            class2=class(obj2);
            if class1=="double" && class2=="FilType"
                r=obj2;
                r.single=operator(obj1,obj2.single);
                if class(obj2.double)=="Filament"
                    r.double.a=operator(obj1,obj2.double.a);
                    r.double.b=operator(obj1,obj2.double.b);
                else
                    r.double=operator(obj1,obj2.double);
                end
                if class(obj2.dimer)=="Filament"
                    r.dimer.a=operator(obj1,obj2.dimer.a);
                    r.dimer.b=operator(obj1,obj2.dimer.b);
                else
                    r.dimer=operator(obj1,obj2.dimer);
                end
                return
            elseif class2=="FilType" && class1=="FilType"
                r=obj2;
                r.single=operator(obj1.single,obj2.single);
                if class(r.double)=="Filament" && class(obj1.double)=="Filament"
                    r.double.a=operator(obj1.double.a,obj2.double.a);
                    r.double.b=operator(obj1.double.b,obj2.double.b);
                else
                    r.double=operator(obj1.double,obj2.double);
                end
                if class(r.dimer)=="Filament" && class(obj1.dimer)=="Filament"
                    r.dimer.a=operator(obj1.dimer.a,obj2.dimer.a);
                    r.dimer.b=operator(obj1.dimer.b,obj2.dimer.b);
                else
                    r.dimer=operator(obj1.dimer,obj2.dimer);
                end
                return
            elseif class2=="double" && class1=="FilType"
                r=obj1;
                r.single=operator(obj1.single,obj2);
                if class(obj1.double)=="Filament"
                    r.double.a=operator(obj1.double.a,obj2);
                    r.double.b=operator(obj1.double.b,obj2);
                else
                    r.double=operator(obj1.double,obj2);
                end
                if class(obj1.dimer)=="Filament"
                    r.dimer.a=operator(obj1.dimer.a,obj2);
                    r.dimer.b=operator(obj1.dimer.b,obj2);
                else
                    r.dimer=operator(obj1.dimer,obj2);
                end
            else
                error("inputs must be of type FilType or double")
            end
        end
    end

    methods(Access={?Formin})
        function obj=add_fils(obj)
            obj.double=obj.double.a+obj.double.b;
            obj.dimer=obj.dimer.a+obj.dimer.b;
        end
    end
end