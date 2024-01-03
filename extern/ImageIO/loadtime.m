function time = loadtime(filename)
%LOAD_TIME This function load the absolute time(world time) from file
% input:
%   - filename: the microscope raw data, *.ims, *.nd2 are supported
% output:
%   - time: the time points relative to the acquisition start, which
%           is the capture time of last plane in a volume
%
%   see also: bfInitLogging, bfGetReader, java.util.Hashtable, datetime

% Copyright (c) 2021, Weihan Li
% SAVETIFF: Version: 1.0.0
%   *** relative time points support
%   *** tif file half supported (bad time loading)
%   *** modified the time output

arguments
    filename (1,1) string = "";
end

time = [];

if ~exist(filename,"file")
    [file,path] = uigetfile({'*.ims','Imaris Files (*.ims)';'*.nd2','Nikon Files (*.nd2)'},...
        'microscope data selector');
    if isequal(file,0)
        warning('microscope data load failed');
        return;
    else
        filename = string(fullfile(path,file));
    end
end

[~,~,ext] = fileparts(filename);
assert(ismember(ext,[".ims",".nd2",".tif"]),"unsupported file format");

% load time info
t = [];
status = bfCheckJavaPath(1);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

% Initialize logging
bfInitLogging();

% Get the channel filler
r = bfGetReader(filename.char(), 0);

% get global meta data
global_meta_data = r.getGlobalMetadata();

% get packaged meta data
store_meta_data = r.getMetadataStore();

% get frames
frames = store_meta_data.getPixelsSizeT(0).getValue();

if ext == ".tif" || ext == ".ims"
    % get delta time relative to world time
    if ~isempty(global_meta_data.get("AcquisitionStart")) ...
            && ~isempty(global_meta_data.get("RecordingDate"))
        dt = seconds(diff([datetime(global_meta_data.get("AcquisitionStart")),...
            datetime(global_meta_data.get("RecordingDate"))]));
    else
        warning('load_time:timeNotFound',"The acquire time not found.");
        return;
    end

    if isa(global_meta_data,'java.util.Hashtable')
        keys = global_meta_data.keys();
        n = 1;
        while keys.hasMoreElements()
            tp_name = ['TimePoint',num2str(n)];
            if global_meta_data.containsKey(tp_name)
                t = [t,datetime(global_meta_data.get(tp_name))]; %#ok<AGROW>
                n = n + 1;
            end
            keys.nextElement();    % inner iterator increase
        end
    else
        error('invalid input format');
    end

    % datetime array t, change it to seconds
    time = [0,cumsum(seconds(diff(t)))]';
    time(frames+1:end) = [];
    time(1) = (time(end)-time(2))/(frames - 2 + eps);   % average estimation
    time = time + dt;

    % control the precision: ms
    time = round(time, 3);
elseif ext == ".nd2"
    slices = store_meta_data.getPixelsSizeZ(0).getValue();
    channels = store_meta_data.getPixelsSizeC(0).getValue();
    time = nan(frames, 1);
    for k = 1:frames-1
        time(k+1) = store_meta_data.getPlaneDeltaT(0,slices*channels*k-1).value();
    end
    time = fillmissing(time, "linear");
    time = round(time, 3);
end

r.close();
end

