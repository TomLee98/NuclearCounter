function volume = bfopen_reg(filename, wbar_flag)
%BFOPEN_REG This function read the image file and load the first series to
%memory
% Input:
%   - filename: the imaging file name
%   - wbar_flag: the flag for viewing processing bar
% Output:
%   - volume: the M*N*P imaging matrix

% Toggle the autoloadBioFormats flag to control automatic loading
% of the Bio-Formats library using the javaaddpath command.
%
% For static loading, you can add the library to MATLAB's class path:
%     1. Type "edit classpath.txt" at the MATLAB prompt.
%     2. Go to the end of the file, and add the path to your JAR file
%        (e.g., C:/Program Files/MATLAB/work/bioformats_package.jar).
%     3. Save the file and restart MATLAB.
%
% There are advantages to using the static approach over javaaddpath:
%     1. If you use bfopen within a loop, it saves on overhead
%        to avoid calling the javaaddpath command repeatedly.
%     2. Calling 'javaaddpath' may erase certain global parameters.
autoloadBioFormats = 1;

% Toggle the stitchFiles flag to control grouping of similarly
% named files into a single dataset based on file numbering.
stitchFiles = 0;

% To work with compressed Evotec Flex, fill in your LuraWave license code.
%lurawaveLicense = 'xxxxxx-xxxxxxx';

% -- Main function - no need to edit anything past this point --

% load the Bio-Formats library into the MATLAB environment
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

% Prompt for a file if not input
if nargin == 0 || exist(filename, 'file') == 0
  [file, path] = uigetfile(bfGetFileExtensions, 'Choose a file to open');
  filename = [path file];
  if isequal(path, 0) || isequal(file, 0), return; end
end

% Initialize logging
bfInitLogging();

% Get the channel filler
r = bfGetReader(filename, stitchFiles);

% Test plane size
planeSize = javaMethod('getPlaneSize', 'loci.formats.FormatTools', r);
if planeSize/(1024)^3 >= 2
    error(['Image plane too large. Only 2GB of data can be extracted '...
        'at one time. You can workaround the problem by opening '...
        'the plane in tiles.']);
end
numSeries = r.getSeriesCount();
assert(numSeries>=1,"series lost");

info = r.getMetadataStore();

width = info.getPixelsSizeX(0).getValue();
height = info.getPixelsSizeY(0).getValue();
pixelType = r.getPixelType();
bpp = javaMethod('getBytesPerPixel', 'loci.formats.FormatTools', ...
    pixelType);
numImages = r.getImageCount();

volume = zeros(height, width, numImages, "uint"+string(8*bpp));

fprintf('Reading series #%d\n', 1);

if wbar_flag == true
    bar = waitbar(0,"Start reading data...");
  
    r.setSeries(0);

    for slice_k = 1:numImages
        volume(:,:,slice_k) = bfGetPlane(r, slice_k);
        switch mod(slice_k, 100)
            case 1
                procs = slice_k/numImages;
                waitbar(procs,bar,['loading process: ',num2str(100*procs,'%.1f'),' %']);
            otherwise
        end
    end

    close(bar);
else
    r.setSeries(0);

    for slice_k = 1:numImages
        volume(:,:,slice_k) = bfGetPlane(r, slice_k);
    end
end

r.close();
end

