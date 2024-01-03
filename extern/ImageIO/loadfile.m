function [opts,rt,movie,filename] = loadfile(filename,omitif,order_out)
%LOADFILE: This function for loading data (.ims,.tif,.nd2)
% input:
%   - filename:char array or string, the movie file name (with full path), can be empty
%   - omitif: the flag for omit tiff file information lost, false as default
%   - order_out: the output stacks order, ["X","Y","C","Z","T"] as default
% output:
%   - opts: some useful movie information, 1 X 12 table, with variables:
%           width   height  channels slices  frames images  xRes    yRes
%           zRes    bitDepth    dimOrder    cOrder
%   - movie: 5-D or 4-D grayscale matrix, type as source
%   - rt: the relative camera time, n*1 array, seconds as unit
%   - filename: the file name, char array
%
%   see also: imload

% Copyright (c) 2022-2023, Weihan Li
% LOADFILE: 
%
% Verison: 1.2.0
%   *** using faster library for *.ims data loading


arguments
    filename string;
    omitif (1,1) logical = false;
    order_out (1,5) string = ["X","Y","C","Z","T"];
end

nargoutchk(1, 4);

opts = [];
movie = [];
rt = [];

% check whether file is existing or not
if isempty(filename)
    % uigetfile seletion
    [file,path] = uigetfile(...
        {'*.ims;*.tif;*.nd2','volume files (*.ims,*.tif,*.nd2)'},...
        'volume data selector');
    if isequal(file,0) || isequal(path,0)
        return;
    else
        filename = fullfile(path,file);
    end
end

if nargout < 3
    [status, info] = imload(filename, true, omitif);
else
    [status, info, img] = imload(filename, true, omitif);
end

if status == -1
    error("Loading failed.");
else
    opts = info.opts;
    rt = info.rt;

    if nargout >= 3
        % using order_out for transformed order
        movie = reorder(img, opts.dimOrder, order_out);
        disp("Loading success.");
    end
end


    function mov = reorder(mov, order_old, order_new)
        [~, order] = ismember(order_new, order_old);
        mov = permute(mov, order);
    end
end
