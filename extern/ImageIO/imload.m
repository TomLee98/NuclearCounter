function [status, info, img] = imload(file, wbar_flag, omitif)
%IMLOAD This function check the file by using the ext and select the
%reader for reading images into memory
% Input:
%   - file: string, the file full path
%   - wbar_flag: boolean, the waitbar showing flag, true as default
% Output:
%   - status: the read status
%   - info: 1*1 struct, with the file info [opts](table) and real time 
%           stamp [rt](vector)
%   - img: the image stack data, 5-D matrix

% Copyright (c) 2022-2024, Weihan Li
% SAVEFILE: 
% Version: 1.0.0
%   *** basic saving function
%   *** status indicate the function status
%   *** many export data format support
% Version: 1.1.0
%   *** import *.ims support library for fast reading
%   *** import *.nd2 support library(extern C) for fast reading
%   *** import *.tiff support library(extern Python) for fast reading

arguments
    file (1,1) string {mustBeFile};
    wbar_flag (1,1) logical = true;
    omitif (1,1) logical = false;
end

[~, ~, ext] = fileparts(file);

if nargout == 3
    switch ext
        case ".ims"
            [status, info, img] = loadImsFile(file, wbar_flag);
        case ".nd2"
            [status, info, img] = loadNd2File(file, wbar_flag);
        case ".tif"
            [status, info, img] = loadTiffFile(file, wbar_flag, omitif);
        otherwise
            error("Unsupported file format.");
    end
elseif nargout == 2
    switch ext
        case ".ims"
            [status, info] = loadImsFile(file, wbar_flag);
        case ".nd2"
            [status, info] = loadNd2File(file, wbar_flag);
        case ".tif"
            [status, info] = loadTiffFile(file, wbar_flag, omitif);
        otherwise
            error("Unsupported file format.");
    end
else
    error("Invalid output arguments number.");
end

end


function [status, info, img] = loadImsFile(file, wbar_flag)
% generate HDF5 reader
GID = H5G.open(H5F.open(file), '/DataSet');
dataobj = DatasetReader(GID);

nargoutchk(2,3);

% extract the information
width = dataobj.SizeX;
height = dataobj.SizeY;
channels = dataobj.SizeC;
slices = dataobj.SizeZ;
frames  = dataobj.SizeT;
images = channels*slices*frames;
xRes = (dataobj.ExtendMaxX-dataobj.ExtendMinX)/dataobj.SizeX;   % um/pix
yRes = (dataobj.ExtendMaxY-dataobj.ExtendMinY)/dataobj.SizeY;   % um/pix
zRes = (dataobj.ExtendMaxZ-dataobj.ExtendMinZ)/(dataobj.SizeZ-1);   % um/pix
dataType = string(dataobj.DataType);
if dataType == "float"
    dataType = "single";
end
dimOrder = ["X","Y","Z","C","T"];
cOrder = getChannelOrder(dataobj.ChannelInfo);

% generate opts
opts = table(width,...          % #pixel
            height,...          % #pixel
            channels,...        % #channels
            slices,...          % #slices
            frames,...          % #volumes
            images,...          % #images
            xRes,...            % μm/pixel
            yRes,...            % μm/pixel
            zRes,...            % μm/layer z
            dataType,...        % uint8, uint16, single
            dimOrder,...        % dimention order array
            cOrder);            % color channel order

% extract time stamps
rt = getRelativeTime(dataobj.Timestamps);

info.opts = opts;
info.rt = rt;

status = 0;

if nargout == 3
    % load the image data
    img = dataobj.GetData(wbar_flag);
end

    function cOrder = getChannelOrder(cInfo)
        order = ["r","g","b"];
        cOrder = strings(1, numel(cInfo));
        for k = 1:numel(cInfo)
            [~, p] = max(cInfo(k).Color);
            cOrder(k) = order(p);
        end
    end

    function rt = getRelativeTime(timeStamp)
        rt = datetime(string(timeStamp),"Format","uuuu-MM-dd HH:mm:ss.SSS");
        rt = seconds(rt - rt(1));
    end
end

function [status, info, img] = loadNd2File(file, wbar_flag)

nargoutchk(2,3);

file = char(file);

% we need to use the bfmatlab package
status = bfCheckJavaPath(1);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);
bfInitLogging();
r = bfGetReader(file, 0);
info_bf = r.getMetadataStore();

% extract the information
width = info_bf.getPixelsSizeX(0).getValue();
height = info_bf.getPixelsSizeY(0).getValue();
channels = info_bf.getPixelsSizeC(0).getValue();
slices = info_bf.getPixelsSizeZ(0).getValue();
frames = info_bf.getPixelsSizeT(0).getValue();
images = channels*slices*frames;
xRes = double(info_bf.getPixelsPhysicalSizeX(0).value);
yRes = double(info_bf.getPixelsPhysicalSizeY(0).value);
if slices ~= 1
    zRes = double(info_bf.getPixelsPhysicalSizeZ(0).value);
else
    zRes = 0;   % only one plane, no resolution at z direction
end
dimOrder = string(split(info_bf.getPixelsDimensionOrder(0).getValue(),""));
dimOrder = dimOrder(strlength(dimOrder)>0)';
cOrder = getChannelOrder(info_bf);
dataType = string(info_bf.getPixelsType(0).getValue());
if dataType == "float"
    dataType = "single";
end

% generate opts
opts = table(width,...          % #pixel
            height,...          % #pixel
            channels,...        % #channels
            slices,...          % #slices
            frames,...          % #volumes
            images,...          % #images
            xRes,...            % μm/pixel
            yRes,...            % μm/pixel
            zRes,...            % μm/layer z
            dataType,...        % uint8, uint16, single
            dimOrder,...        % dimention order array
            cOrder);            % color channel order

rt = getRelativeTime(r);

r.close();

info.opts = opts;
info.rt = rt;

status = 0;

if nargout == 3
    if ispc()
        % only PC support fast loading
        img = nd2Open3DVolume_reg(file, wbar_flag);
    elseif isunix()
        img = bfOpen3DVolume_reg(file, wbar_flag);
    end

     % reconstruct the image stack
    img = imreshape(img, opts);
end

    function mov = nd2Open3DVolume_reg(file, wbar_flag)
        if wbar_flag == true
            fig = uifigure("Visible","off","WindowStyle","modal",...
                "Resize","off");
            fig.Position(3:4) = [300, 75];
            set(fig, "Visible", "on");
            uiprogressdlg(fig,'Indeterminate','on', ...
                'Message', '        loading...','Icon','info',...
                "Interpreter","tex");
        end

        % use nd2 library load total volumes
        [FilePointer, ImagePointer, ImageReadOut] = ND2Open(file);
        numImages = calllib('Nd2ReadSdk', 'Lim_FileGetSeqCount', FilePointer);
        Num=1:numImages;
        % mov is a c-by-1 cell, c for channels number
        mov = ND2Read(FilePointer, ImagePointer, ImageReadOut, Num);
        ND2Close(FilePointer)
        clear("FilePointer","ImagePointer","ImageReadOut");

        for n = 2:numel(mov)
            mov{1} = cat(4, mov{1}, mov{n});
            mov{n} = [];    % free memory
        end
        mov = mov{1};
        % permute as (Y,X,C,Z(T)): nd2 file fixed dimension order
        mov = permute(mov, [1,2,4,3]);

        if wbar_flag == true
            delete(fig);
        end
    end

    function cOrder = getChannelOrder(info)
        cc = info.getChannelCount(0);
        C = zeros(cc,4,'uint8');
        for k = 0:cc-1
            cinfo = info.getChannelColor(0,k);
            % some format there is no color channel information such as
            % tiff file exported from fiji
            if isempty(cinfo)
                cOrder = "";
                return;
            end
            C(k+1,:) = [cinfo.getRed(),cinfo.getGreen(),cinfo.getBlue(),...
                cinfo.getAlpha()];
        end
        C_na = C(:,1:3);  % cut the color map without alpha channel
        % select the maximum value to represent the channel
        [~,order] = max(C_na,[],2);
        cOrder = ["r","g","b"];
        cOrder = cOrder(order);
    end

    function rt = getRelativeTime(bf_reader)
        % get packaged meta data
        store_meta_data = bf_reader.getMetadataStore();
        rt = nan(frames, 1);
        for k = 1:frames-1
            rt(k+1) = store_meta_data.getPlaneDeltaT(0,slices*channels*k-1).value();
        end
        rt = fillmissing(rt, "linear");
        rt = round(rt, 3);
    end

    function mov = imreshape(mov, opts)
        bindingArray = ["X",    "Y",    "C",    "Z",    "T";...
            "height","width","channels","slices","frames"];
        [~,Loc] = ismember(opts.dimOrder,bindingArray(1,:));
        mov = reshape(mov,...
            opts.(bindingArray(2,Loc(1))),...
            opts.(bindingArray(2,Loc(2))),...
            opts.(bindingArray(2,Loc(3))),...
            opts.(bindingArray(2,Loc(4))),...
            opts.(bindingArray(2,Loc(5))));
    end
end

function [status, info, img] = loadTiffFile(file, wbar_flag, omitif)
nargoutchk(2,3);

file = char(file);

% we need to use the bfmatlab package
status = bfCheckJavaPath(1);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);
bfInitLogging();
r = bfGetReader(file, 0);
info_bf = r.getMetadataStore();

% extract the information
width = info_bf.getPixelsSizeX(0).getValue();
height = info_bf.getPixelsSizeY(0).getValue();
channels = info_bf.getPixelsSizeC(0).getValue();
slices = info_bf.getPixelsSizeZ(0).getValue();
frames = info_bf.getPixelsSizeT(0).getValue();
images = channels*slices*frames;
dimOrder = string(split(info_bf.getPixelsDimensionOrder(0).getValue(),""));
dimOrder = dimOrder(strlength(dimOrder)>0)';
cOrder = getChannelOrder(info_bf);
dataType = string(info_bf.getPixelsType(0).getValue());
if dataType == "float"
    dataType = "single";
end

% ask for resolution information
if omitif == false
    % Ask the user to input necessary imformation
    prompt = {'Enter the z scan thickness(\mum):',...
        'Enter the binning size:',...
        'Enter total optical magnification:',...
        'Enter the camera pixel size(\mum):'};
    dlgtitle = 'Microscope parameters setting';
    dims = [1 40];
    definput = {'1.5';'2';'60';'6.5'};
    dlgopts.Interpreter = 'tex';
    dlgout = inputdlg(prompt,dlgtitle,dims,definput,dlgopts);
    if isempty(dlgout), dlgout = definput; end

    % extract and calculate the X,Y resolution
    zRes = str2double(dlgout{1});
    binsz = str2double(dlgout{2});
    tcom = str2double(dlgout{3});
    cps = str2double(dlgout{4});

    % since the camera physical pixel is square,
    % and the binning size is also square
    xRes = cps/tcom*binsz;
    yRes = xRes;
else
    % the resolution information is useless
    xRes = 1.0;
    yRes = 1.0;
    zRes = 1.0;
end


% generate opts
opts = table(width,...          % #pixel
            height,...          % #pixel
            channels,...        % #channels
            slices,...          % #slices
            frames,...          % #volumes
            images,...          % #images
            xRes,...            % μm/pixel
            yRes,...            % μm/pixel
            zRes,...            % μm/layer z
            dataType,...        % uint8, uint16, single
            dimOrder,...        % dimention order array
            cOrder);            % color channel order

rt = getRelativeTime(r);

r.close();

info.opts = opts;
info.rt = rt;

status = 0;

if nargout == 3
    % check the environment
    pyflag = isPyReady();

    if pyflag == true
        img = pyOpen3DVolume_reg(file, wbar_flag);
    else
        img = bfOpen3DVolume_reg(file, wbar_flag);
    end

    % reconstruct the image stack
    img = imreshape(img, opts);
end

    function mov = pyOpen3DVolume_reg(file, wbar_flag)
        % check the python environment path
        if count(py.sys.path,'/extern/ImageIO/TIFFRW/load_tiff.py') == 0
            insert(py.sys.path,int32(0), ...
                '/extern/ImageIO/TIFFRW/load_tiff.py');
        end
        fname = py.str(file);
        if wbar_flag == true
            fig = uifigure("Visible","off","WindowStyle","modal",...
                "Resize","off");
            fig.Position(3:4) = [300, 75];
            set(fig, "Visible", "on");
            uiprogressdlg(fig,'Indeterminate','on', ...
                'Message', '        loading...','Icon','info',...
                "Interpreter","tex");
        end
        try
            % imagej stack: 'TZCYXS'
            mov = pyrunfile("load_tiff.py", "vol", file=fname);
        catch ME
            throwAsCaller(ME);
        end
        % convert img from ndarray to matlab value
        mov = cast(mov, dataType);
        if ndims(mov) == 3
            mov = permute(mov, [2,3,1]);   % to (Y,X,Z)
        elseif ndims(mov) == 4
            mov = permute(mov, [3,4,2,1]); % to (Y,X,Z,T)
        elseif ndims(mov) == 5
            mov = permute(mov, [4,5,3,2,1]); % to (Y,X,C,Z,T)
        else
            throw(MException("imload:invalidImagesStackDimension", ...
                "Image stack dimension > 5 is not supported."));
        end

        if wbar_flag == true && isvalid(fig)
            delete(fig);
        end
    end

    function cOrder = getChannelOrder(info)
        cc = info.getChannelCount(0);
        C = zeros(cc,4,'uint8');
        for k = 0:cc-1
            cinfo = info.getChannelColor(0,k);
            % some format there is no color channel information such as
            % tiff file exported from fiji
            if isempty(cinfo)
                disp("The channel order not found. Manual inputs are required.");
                cn = info.getPixelsSizeC(0).getValue(); % color channels number
                cdlgtitle = 'Camera color channel setting';
                cdims = [1 40];
                switch cn
                    case 1
                        cprompt = {'Enter the channel color(r/g/b):'};
                        cdefinput = {'g'};
                    case 2
                        cprompt = {'Enter the channel-1 color(r/g/b):',...
                                   'Enter the channel-2 color(r/g/b):'};
                        cdefinput = {'r','g'};
                    otherwise
                        error("Unsupported Color Channel.");
                end
                cdlgout = inputdlg(cprompt,cdlgtitle,cdims,cdefinput);
                if isempty(cdlgout) || numel(cdlgout) ~= cn
                    cdlgout = cdefinput; 
                end
                cOrder = reshape(string(cdlgout),1,[]);

                return;
            end
            C(k+1,:) = [cinfo.getRed(),cinfo.getGreen(),cinfo.getBlue(),...
                cinfo.getAlpha()];
        end
        C_na = C(:,1:3);  % cut the color map without alpha channel
        % select the maximum value to represent the channel
        [~,order] = max(C_na,[],2);
        cOrder = ["r","g","b"];
        cOrder = cOrder(order);
    end

    function rt = getRelativeTime(bf_reader)
        % get packaged meta data
        global_meta_data = bf_reader.getGlobalMetadata();

        % get delta time relative to world time
        if ~isempty(global_meta_data.get("AcquisitionStart")) ...
                && ~isempty(global_meta_data.get("RecordingDate"))
            dt = seconds(diff([datetime(global_meta_data.get("AcquisitionStart")),...
                datetime(global_meta_data.get("RecordingDate"))]));
        else
            warning("The acquire time not found. Frame Indices replaced.");
            rt = [];
            return;
        end

        if isa(global_meta_data,'java.util.Hashtable')
            keys = global_meta_data.keys();
            n = 1;
            t = [];
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
        rt = [0,cumsum(seconds(diff(t)))]';
        rt(frames+1:end) = [];
        if frames > 1
            rt(1) = (rt(end)-rt(2))/(frames - 2 + eps);   % average estimation
            rt = rt + dt;
        end

        % control the precision: ms
        rt = round(rt, 3);
    end

    function mov = imreshape(mov, opts)
        bindingArray = ["X",    "Y",    "C",    "Z",    "T";...
                        "height","width","channels","slices","frames"];
        [~,Loc] = ismember(opts.dimOrder,bindingArray(1,:));
        mov = reshape(mov,...
            opts.(bindingArray(2,Loc(1))),...
            opts.(bindingArray(2,Loc(2))),...
            opts.(bindingArray(2,Loc(3))),...
            opts.(bindingArray(2,Loc(4))),...
            opts.(bindingArray(2,Loc(5))));
    end
end
