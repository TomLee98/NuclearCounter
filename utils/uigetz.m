function roi = uigetz(range)
%UIGETZ This function get z slices and return array
% input:
%   - range: 1-by-2 double, the z limits
% output:
%   - roi: n-by-1 positive integer array

arguments
    range (1,2) double {mustBeInteger, mustBePositive}
end

assert(range(1) < range(2), "uigetz:invalidZRange", ...
    "First location must be less than second location.");

prompt = {'Enter z range(split with comma):'};
dlgtitle = "crop z";
fieldsize = [1 45];
definput = {'1:end'};

answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

if ~isempty(answer)
    answer{1} = string(answer{1}).replace("end", string(range(2))).char();
    roi = unique(str2num(answer{1}),"rows","sorted"); %#ok<ST2NM>
    assert((roi(1)>=range(1))&&(roi(end)<=range(2)&&(roi(1)<roi(end))), ...
        "uigetz:invalidZRange", "Location indices must be in range.");
else
    roi = [];
end
end

