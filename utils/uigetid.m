function id = uigetid(id, tl)
%UIGETID This function uses inputdlg to get identity by user input
% input:
%   - id: the initial id, double, positive integer
%   - tl: the input dialog title
% output:
%   - id: the output id

arguments
    id (1,1) double {mustBePositive, mustBeInteger} = 1;
    tl (1,1) string  = "rename";
end

prompt = {'Enter new label:'};
dlgtitle = tl;
fieldsize = [1 45];
definput = {num2str(id)};

answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
if ~isempty(answer)
    id = str2double(answer{1});
else
    id = [];
end
end

