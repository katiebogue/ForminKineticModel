function output = runpyfunction(options,dependencies,script,func,inputs)
% RUNPYFUNCTION run a python function from a user defined python script in
% the specified directory
    %   output=
    %   RUNPYFUNCTION(options,dependencies,script,func,inputs) runs the
    %   python function func from script with the specified inputs
    %
    %   Inputs:
    %         options      : structure/class containing the location of the 
    %                        python script in the element "python_path"
    %         dependencies : additonal python dependencies to load (must be
    %                        a cell array)
    %         script       : the python script containing the target 
    %                        function
    %         func         : the python function to run
    %         inputs       : cell array of inputs to the python function
    %   
    %   Outputs:
    %         output : the output of the python function (may be a python
    %         tuple)
    %   
    %   Example:
    %       to run get_formin_info.py:
    %           runpyfunction(opt,{'bioservices'},'get_formin_info','main',{"-1", 'sequence', 4, 1, 1, 1, 5, 3, 1})
    %
    %   See also DOUBLE, CELL.
py.sys.path().append(options.python_path);
for dep=1:length(dependencies)
    py.importlib.import_module(dependencies{dep});
end
pyscript=py.importlib.import_module(script);
output= pyscript.(func)(inputs{:});

end