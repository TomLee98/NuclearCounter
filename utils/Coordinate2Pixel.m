function [R_pix,W_pix,H_pix] = Coordinate2Pixel(R,rdx,rdy)
% COORDINATE2PIXEL This function transform a rectangle in the nature coordinate 
% system to pixel coordinate system with linear mapping
% input:
%   - R: the rectangle matrix in nature coordinate system [x1,y1;x2,y2]
%   - rdx: metric unit per pixel at width direction
%   - rdy: metric unit per pixel at height direction
% output:
%   - R_pix: the rectangle matrix in pixel coordinate system [x1,y1;x2,y2]
%   - H_pix: the height of rectangle in pixel coordinate system
%   - W_pix: the width of rectangle in pixel coordinate system

% Version 1.1.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    R (2,2) double;     % [x1,y1;x2,y2]
    rdx (1,1) double;   % pix^-1, dx per pixel
    rdy (1,1) double;   % pix^-1, dy per pixel
end

% convert coordinate to pixel size
W_pix = ceil(R(2,1)/rdx)-fix(R(1,1)/rdx);
H_pix = ceil(R(2,2)/rdy)-fix(R(1,2)/rdy);
R_pix = [fix(R(1,1)/rdx)+1,fix(R(1,2)/rdy)+1;...
    ceil(R(2,1)/rdx),ceil(R(2,2)/rdy)];
end

