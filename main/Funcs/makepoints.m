function [colors,shapes]= makepoints(len,background_color)
% MAKEPOINTS  Generates array of distinguishable colors and (repeating)
%             shapes. 
    %
    %   C = MAKEPOINTS(N,color) outputs array of N colors distinguishable 
    %   from color and themselves and array of 8 shapes repeated ceil(N/8)+1
    %   times.
    %   
    %   Inputs:
    %       len              : desired length of output color & shape array 
    %       background_color : color for all colors to be distinct from
    %
    %   Outputs:
    %       1) colors : array of distinguishable colors of input length  
    %       2) shapes : array of 8 shapes repeated ceil(len/8)+1 times
    %   
    %   See also OPTIONS, FIGURESAVE.

    %create colors matrix
    rgbmatrix = distinguishable_colors(len,background_color);
    colors = rgb2hex(rgbmatrix);
    colors=string(colors);
    
    %create matrix of plot shapes
    shapes = ['o';'s';'d';'^';'v';'<';'p';'h'];
    shapes = repmat(shapes,ceil(len/8)+1,1);
end