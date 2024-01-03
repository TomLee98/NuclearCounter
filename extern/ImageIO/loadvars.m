function [opts,rt,movie,varnames] = loadvars()
%LOADVARS This function call VarSelector for loading variables in base
% workspace
% input:
% 
% output:
%   - opts: some useful movie information, 1 X 12 table, with variables:
%           width   height  channels slices  frames images  xRes    yRes
%           zRes    bitDepth    dimOrder    cOrder
%   - movie: 5-D or 4-D grayscale matrix, type as source
%   - rt: the relative camera time, n*1 array, seconds as unit
%   - varnames: the selected var names, cell

opts = [];
movie = [];
rt = [];
varnames = [];

vs_handle = VarSelecter_GUI();

waitfor(vs_handle,"is_closing",true);

% take the information and delete app
info = vs_handle.UIFigure.UserData;
delete(vs_handle);

s = validateinfo(info);
if s ~= 0
    % bad selecting, return with empty array
    return;
end

opts = evalin("base", info.opts);
validateopts(opts);

if opts.channels == 1
    % each movie format as (X,Y,Z,T)
    mov_r = evalin("base", info.mov_r);
    validatemovie(mov_r, 1);
    mov_g = evalin("base", info.mov_g);
    validatemovie(mov_g, 1);
    movie = cat(5,mov_r,mov_g);   % format order as R, G, B
    clearvars mov_r mov_g;
    % (X,Y,Z,T,C) -> (X,Y,C,Z,T)
    movie = permute(movie,[1,2,5,3,4]);

    % modified opts.channels, opts.images, opts.dimOrder, opts.cOrder
    opts.channels = 2;
    opts.images = opts.images * 2;
    opts.dimOrder = ["X","Y","C","Z","T"];
    opts.cOrder = ["r","g"];
else
    movie = evalin("base", info.mov);
    validatemovie(movie, 2);
end

rt = evalin("base", info.rt);
validatert(rt);

varnames = info.varnames;

end

function status = validateinfo(info)
% This function validates the data information
status = 0;

assert(isstruct(info),"validateinfo:InfoFormatError");

f = fieldnames(info);
valid_field = ((numel(f)==4) && all(ismember(f,{'mov','opts','rt','varnames'})))...
    || ((numel(f)==5) && all(ismember(f,{'mov_r','mov_g','opts','rt','varnames'})));

assert(valid_field,'validateinfo:InfoFieldError');

if any(info.varnames == "none")
    status = -1;
end

end

function validateopts(opts)
% This function validates the movie options

if isempty(opts)
    return;
end

assert(istable(opts),"validateopts:OptsIsNotTable");

assert(isequal(size(opts),[1,12]),"validateopts:InvalidSizeofOpts");

v = opts.Properties.VariableNames;
valid_field = all(ismember(v,{'width','height','channels','slices','frames',...
    'images','xRes','yRes','zRes','dataType','dimOrder','cOrder'}));

assert(valid_field,"validateopts:InvalidOptsVariable");

end

function validatemovie(movie, c)
% This function validates the movie options

if isempty(movie)
    return;
end

assert(isnumeric(movie),"validatemovie:InvalidType");

if c == 1
    assert(ndims(movie)==4,"validatemovie:InvalidDimension");
elseif c == 2
    assert(ndims(movie)==5,"validatemovie:InvalidDimension");
    % for checking if dimension order is (-,-,C,-,-)
    assert(size(movie,3)==2,"validatemovie:InvalidFormat");
else
    error("validatemovie:BadCaller","Invalid range of input argument c.");
end

end

function validatert(rt)
% This function validates the real time (camera time)

assert(isnumeric(rt),"validatert:InvalidType");

assert(isvector(rt),"validatert:InvalidFormat");

assert(issorted(rt,"ascend"),"validatert:InvalidOrder");

end
