function status = imsave(file, data, info)
%IMSAVE This function write data to disk. Note that if there is python
%environment, calling tifffile can be faster than OME/bioformat
% Input:
%   - file: string, the file full path
%   - data: 5-D matrix
%   - info: the file info, struct
% Output:
%   - status: the write status

% Copyright (c) 2022-2024, Weihan Li
% IMSAVE: 
% Version: 1.0.0
%   *** basic saving function


[~, ~, ext] = fileparts(file);

switch ext
    case ".tif"
    case ".png"
    case ".gif"
end

end

function saveTiffFile()

end

function savePngFile()

end

function saveGifFile()

end