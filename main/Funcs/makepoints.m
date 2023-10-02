function [colors,shapes]= makepoints(len,background_color)
% MAKEPOINTS  Generates array of distinguishable colors and (repeating)
%             shapes. 
    %
    %   C = MAKEPOINTS(N,color) outputs array of N colors distinguishable from
    %   color and themselves and array of 8 shapes repeated N+1 times.
    %   
    %   Inputs:
    %       len              : desired length of output color array 
    %       background_color : color for all colors to be distinct from
    %
    %   Outputs:
    %       1) colors : array of distinguishable colors of input length  
    %       2) shapes : array of 8 shapes repeated len+1 times
    %   
    %   See also .

    %create colors matrix
    rgbmatrix = distinguishable_colors(len,background_color);
    colors = [];
    for k = 1:length(rgbmatrix(:,1))
        hexcode = rgb2hex(rgbmatrix(k,:));
        hexcode = convertCharsToStrings(hexcode);
        colors = [colors, hexcode];
    end
    
    %create matrix of plot shapes
    shapes = ['o';'s';'d';'^';'v';'<';'p';'h'];
    for k = 1:len
        shapes1 = ['o';'s';'d';'^';'v';'<';'p';'h'];
        shapes = [shapes; shapes1];
    end 
end