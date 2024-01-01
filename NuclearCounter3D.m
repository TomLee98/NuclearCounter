function stat = NuclearCounter3D(filename, varargin)
%NUCLEARCOUNTER3D This function help to count KCs number in a volume
% WARNING: This function require the max resolution image from 1-p confocol
% microscope for best processing: 2048*2048,dz=1\mum
%
%   NuclearCounter3D()
%   stat = NuclearCounter3D(varargin)
%
% input:
%   - filename: the volume file name or [], supported image format: 
%     *.ims, *.nd2, *.tif, note that tiff file need imgopt
%   - (optional) imopts(table): the options of image, usually comes from
%     loadfile, if you pass this argument, inner options will be covered
%   - (optional) bkg: the background remove method, can be "auto" or
%     nonnegtive scalar, 0 for no background need to be removed, 0 as default
%   - (optional) animal(struct):
%       - marker: string scalar or char vector, "mCherry" as default
%       - driver: string scalar or char vector, "OK107" as default
%       - stage: string scalar or char vector, "L1" as default
%       note that write case is omitted
%   - (parameter) nuclear(struct):
%       - radius: the estimation radius of nuclear(gaussian kernel) (\mum),
%       - sigma: the standard deviation of nuclear(g.k.) (\mum)
%       - color: the marker color of nuclear, which could be "r","g","b"
%   - (parameter) imorph(struct):
%       - rehist_quantile: the quantile value for image rehistogram, higher
%         value will omit the darker object
%       - cir_ssy: the sensitivity for imfindcircles, higher may cause more
%         false positive samples
%       - simr_th: the maxium similarity ratio threshold between two circle
%         A and B, simr := #(A & B)/min(#A, #B), 0.8 as default, higher may
%         cause over-segmentation
%       - avgI_th: the average intensity threshold of nuclear, which 
%         should be in range (0, 1), 0.02 as default, higher may cause
%         more false negtive 
%       - gamma: the gamma value in range (0,1) when applying power 
%         transformation for contrast enhancement on dark objects 
%   - (parameter) output(struct):
%       - is_display: flag for displaying or not, true as default
%       - is_saving: flag for saving or not, false as default
%       - format: the image format, "jpg", "tif" or "png", "png" as default
%       - resolution: the graphics resolution(ppi), 150, 300 or 600,
%         300 as default
%   - (parameter) optimal(struct)
%       - max_iter_num: the maximum number of iteration, positive integer,
%         1e3 as default
%       - tol: the tolerance threshold between previous and current loss
%         function, 1e-6 as stop default
%       - state: the optimal switch state, "on" or "off", "off" as default
%       - srange: 1-by-2 double, the nuclears number search range
%       - rlx: the relaxation of symmetric potential well, 2,4,6,...even
%         number
%       - trdof: 1-by-1 double, tradeoff for object finding, from 0 to 1, 
%         0 for more valid objects is important, 1 for less uncertain objects
%
% output:
%   - stat(struct):
%       - imorph(struct) as input but may updated when optimize circles
%         automatically
%       - omorph(struct):
%           - center: (x,y,z) of each nuclear
%           - radius: r of each nuclear
%       - n: number of objects at each plane
%
% see also: imfinfo, bfopen, imbinarize, bwareaopen, imfill, imopen,
% imfindcircles, viscircles, pdist, text, loadfile, fmincon

% Copyright (c) 2022, Weihan Li
% MORPHKC: Version: 1.0.0
% MORPHKC: Version: 1.1.0
%   *** Update the file loading method
%   *** remove single frame loading
% MORPHKC: Version: 1.1.1
%   *** Fix bugs: the compatibility problem with loadfile.m
% MORPHKC: Version: 1.2.0
%   *** continuous KC number replace incremental KC number
%   *** new algorithm, higher accuracy than others
% MORPHKC: Version: 1.3.0
%   *** Keep relative-cross algorithm only, clean up codes
%   *** add parameters auto optimal options
%   *** add image options for batch easily
%   *** add background removed method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
p.StructExpand = false;

default_imopts = [];

default_bkg = 0;

default_animal = struct('marker', "mCherry", ...
                        'driver', "OK107", ...
                        'stage', "L1");

default_imorph = struct('rehist_quantile', 0.975, ...
                        'cir_ssy',0.935, ...
                        'simr_th',0.7, ...
                        'avgI_th',0.015, ...
                        'gamma', 0.55);

default_optimal = struct('maxitn', 200, ...     % for maximum iteration number
                         'tol', 1e-4, ...       % the tolerance
                         'state', 'off', ...    % optimal running flag, 'on' or 'off'
                         'srange', [20,30], ... % the valid search range     
                         'rlx', 2, ...          % the relaxation of symmetric potential well, positive even number
                         'trdof', 0.5);         % the objects tradeoff, 0 ~ 1, 
                                                % 0 for more valid objects is important, 
                                                % 1 for less uncertain objects is important
default_output = struct('is_display', true,...
                        'is_saving', false, ...
                        'validobj_only', false,...
                        'format',   'png', ...
                        'resolution',300);

addRequired(p,'filename', @(x)validateattributes(x, "string", "scalartext"));
addOptional(p, 'imopts', default_imopts);
addOptional(p, 'bkg', default_bkg, @(x)validateattributes(x, {'string', 'numeric'},"scalar"));
addOptional(p,'animal',default_animal);
addParameter(p,'imorph',default_imorph);
addParameter(p,'optimal',default_optimal);
addParameter(p,'output',default_output);
parse(p, filename, varargin{:});

filename = p.Results.filename;
background = p.Results.bkg;
animal = p.Results.animal;
imorph = p.Results.imorph;
output = p.Results.output;
optimal = p.Results.optimal;
imopts = p.Results.imopts;
[optlb, optub, nuclear] = genoptbound(animal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROGRAM PROCESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%
stat = struct('omorph', struct('center',[], 'radius',[]), ...
              'n', [], ...
              'imorph', imorph);

if isempty(imopts)
    [imopts, ~, vol] = loadfile(filename, false);
else
    [~, ~, vol] = loadfile(filename, true);
end
if isempty(vol)
    return;
end

if numel(size(vol)) == 5 || numel(size(vol)) == 2
    error("MORPHKC: Unsupported data dimension.");
else
    if numel(size(vol)) == 4 && size(vol, 3) > 1
        if imopts.frames > 1
            sltframe = nan;
            while isnan(sltframe) || sltframe<1 || sltframe>imopts.frames
                sltframe = input(sprintf("Please input the frame index selected[1-%d]:", ...
                    imopts.frames));
            end
            vol = vol(:,:,nuclear.color==imopts.cOrder,sltframe);
        else
            vol = vol(:,:,nuclear.color==imopts.cOrder,:);
        end
    end
    vol = squeeze(vol);
end

tradeoff_coeff = optimal.trdof;
search_bound = optimal.srange;
rlx_coeff = optimal.rlx;

[stat.omorph.center, stat.omorph.radius, stat.n, v_adjusted, imorph] ...
    = counts(vol,nuclear,imopts,imorph, optimal);

stat.imorph = imorph;

counter_output(v_adjusted, stat.omorph.center, stat.omorph.radius, stat.n, ...
                nuclear, imopts, ...
                'output', output);

    function [center,radius,counts,volumeb,imorph] = counts(vol,nuclear,imopts,imorph, optimal)
        % remove background at first
        if isnumeric(background)
            vol = vol - background;
        elseif isstring(background) || ischar(background)
            switch string(background)
                case "auto"
                    % extract the angle border and calculate the median as the
                    % estimation of background
                    bkvol = vol([1:3, end-2:end], [1:3, end-2:end], :);
                    vol = vol - median(bkvol, "all");
                otherwise
                    error("NuclearCounter3D:invalidBackgroundAlgorithm", ...
                        "Invalid optional argument: bkg");
            end
        end

        % rehist for enhancing contrast
        vol = rehist(vol, imopts, imorph);

        switch optimal.state
            case 'off'
                % use the given parameters run the counter
                [volumeb,center,radius,counts] = ...
                    relative_cross(vol,nuclear,imopts,imorph, true);
            case 'on'
                % optimize the imaging process parameters
                opf = @(x)opfunc(x, vol, nuclear, imopts);

                options = optimoptions("ga",...
                    "MaxGenerations", optimal.maxitn, ...
                    "FunctionTolerance",optimal.tol, ...
                    "ConstraintTolerance", sqrt(optimal.tol), ...
                    "PopulationSize",100, ...
                    "UseParallel", true, ...
                    "MaxTime", 30);

                % using genetic algorithm to solve the global optimal
                % problem
                if isempty(gcp("nocreate"))
                    parpool("Processes");
                end
                
                [mvx, fval] = ga(opf, numel(optlb), [],[],[],[], ...
                    optlb,optub, [], options);

                % delete(gcp("nocreate"));

                % note that gamma as a hidden variable
                imorph = struct('rehist_quantile', mvx(1), ...
                                'cir_ssy',         mvx(2), ...
                                'simr_th',         mvx(3), ...
                                'avgI_th',         mvx(4), ...
                                'gamma',           mvx(5));
                fprintf("Optimal Value: %.2f\n", fval);
                fprintf("Optimized Parameters:\n " + ...
                    "*[Rehist Quantile]     %.3f\n " + ...
                    "*[Circle Sensitivity]  %.3f\n " + ...
                    "*[Merge Threshold]     %.3f\n " + ...
                    "*[Quality Threshold]   %.3f\n ", ...
                    mvx(1), mvx(2), mvx(3), mvx(4));

                % rerun the algorithm
                [volumeb,center,radius,counts] = ...
                    relative_cross(vol,nuclear,imopts,imorph);
        end

        nan_num = sum(cellfun(@(x)sum(isnan(x)), counts));
        subset_m = cellfun(@(x)max(x,[],"all","omitmissing"),counts, ...
            "UniformOutput",false);
        if ~all(cellfun(@(x)isempty(x), subset_m))
            valid_num = max(cell2mat(subset_m));
        else
            valid_num = 0;
        end

        fprintf("Nuclear Objects Found: %d\n " + ...
            "Circles without Belonging: %d\n", ...
            valid_num, nan_num);
    end

    % oci for relative-cross algorithm
    function [c,r,id] = oci(v, vol, r_range, nuclear, imopts, imorph, bar)
        % each c element: [x, y, warn_flag]
        c = cell(size(vol, 3), 1);
        r = cell(size(vol, 3), 1);
        id = cell(size(vol, 3), 1);

        % number of new nulears are found at plane k
        nk = zeros(size(vol, 3), 1);

        % this function using "human-like" method
        vol = padarray(vol,[0,0,1],0,"both");

        [centers_pre, radii_pre] = imfindcircles(vol(:,:,1),r_range,'Sensitivity',imorph.cir_ssy);
        [centers_cur, radii_cur] = imfindcircles(vol(:,:,2),r_range,'Sensitivity',imorph.cir_ssy);

        s = 2;  % slice iterator
        while true
            if s == size(vol, 3)
                break;
            end

            id{s-1} = zeros(size(radii_cur));

            [centers_next, radii_next] = ...
                imfindcircles(vol(:,:,s+1),r_range,'Sensitivity',imorph.cir_ssy);

            scan_code = ~[isempty(centers_pre), isempty(centers_cur), isempty(centers_next)];

            if ismember(scan_code, [false,false,false; ...
                    false,false,true; ...
                    true,false,false; ...
                    true,false,true], "rows")
                % do nothing
                % which means no nuclear was found on current plane
            elseif ismember(scan_code, [false,true,true; ...
                    true,true,true; ...
                    true,true,false], "rows")
                if scan_code(1) == false
                    codes = zeros(size(centers_cur, 1), 4);

                    % calculate the paired distance
                    dist_cur_next = pdist2(centers_cur, centers_next);
                    [d_cnm, idx_cnm] = min(dist_cur_next,[],2);
                    for row = 1:size(codes, 1)
                        % there must be codes(:,1) = 0

                        % find the next plane circle
                        if d_cnm(row) < radii_cur(row) ...
                                && simr(radii_cur(row), radii_next(idx_cnm(row)), d_cnm(row)) >= imorph.simr_th
                            codes(row, 2) = 1;
                        else
                            codes(row, 2) = 0;
                        end

                        % there must be codes(:,3) = 1
                        codes(row, 3) = 1;

                        % calculate the relative size
                        if codes(row, 2) == 1 && ...
                                radii_next(idx_cnm(row)) >= radii_cur(row) - 0.5*nuclear.sigma/imopts.xRes
                            % Note that 95% range for significant
                            % changing detection, split 'decrease no significance'
                            % as 'increase part'
                            codes(row, 4) = 1;
                        else
                            codes(row, 4) = -1;
                        end
                    end
                elseif scan_code(3) == false
                    codes = zeros(size(centers_cur, 1), 4);

                    % calculate the paired distance
                    dist_cur_pre = pdist2(centers_cur, centers_pre);
                    [d_cpm, idx_cpm] = min(dist_cur_pre,[],2);
                    for row = 1:size(codes, 1)
                        % find the previous plane circle
                        if d_cpm(row) < radii_cur(row) ...
                                && simr(radii_cur(row), radii_pre(idx_cpm(row)), d_cpm(row)) >= imorph.simr_th
                            codes(row, 1) = 1;
                        else
                            codes(row, 1) = 0;
                        end
                        % there must be codes(:,2) = 0
                        % calculate the relative size
                        if codes(row, 1) == 0 || ...
                                radii_cur(row) >= radii_pre(idx_cpm(row)) - 0.5*nuclear.sigma/imopts.xRes
                            % Note that 95% range for significant
                            % changing detection, split 'decrease no significance'
                            % as 'increase part'
                            codes(row, 3) = 1;
                        else
                            codes(row, 3) = -1;
                        end
                        % there must be codes(:,4) = -1
                        codes(row, 4) = -1;
                    end
                else
                    % current centers vs. pre and post
                    % generate the codes for center
                    codes = zeros(size(centers_cur, 1), 4);

                    % calculate the paired distance
                    dist_cur_pre = pdist2(centers_cur, centers_pre);
                    dist_cur_next = pdist2(centers_cur, centers_next);

                    [d_cpm, idx_cpm] = min(dist_cur_pre,[],2);
                    [d_cnm, idx_cnm] = min(dist_cur_next,[],2);

                    for row = 1:size(codes, 1)
                        % find the previous plane circle
                        if d_cpm(row) < radii_cur(row) ...
                                && simr(radii_cur(row), radii_pre(idx_cpm(row)), d_cpm(row)) >= imorph.simr_th ...
                                && avg_intensity(s-2, centers_pre(idx_cpm(row),:), radii_pre(idx_cpm(row)), v) > 255*imorph.avgI_th
                            codes(row, 1) = 1;
                        else
                            codes(row, 1) = 0;
                        end
                        % find the next plane circle
                        if d_cnm(row) < radii_cur(row) ...
                                && simr(radii_cur(row), radii_next(idx_cnm(row)), d_cnm(row)) >= imorph.simr_th ...
                                && avg_intensity(s, centers_next(idx_cnm(row),:), radii_next(idx_cnm(row)), v) > 255*imorph.avgI_th
                            codes(row, 2) = 1;
                        else
                            codes(row, 2) = 0;
                        end
                        % calculate the relative size
                        if codes(row, 1) == 0 || ...
                                radii_cur(row) >= radii_pre(idx_cpm(row)) - 0.5*nuclear.sigma/imopts.xRes
                            codes(row, 3) = 1;
                        else
                            codes(row, 3) = -1;
                        end
                        if codes(row, 2) == 1 && ...
                                radii_next(idx_cnm(row)) >= radii_cur(row) - 0.5*nuclear.sigma/imopts.xRes
                            codes(row, 4) = 1;
                        else
                            codes(row, 4) = -1;
                        end
                    end
                end

                % distribute the labels
                circle_status = is_new(codes);
                new_circles_tot = (circle_status == 1);
                new_circles_valid = new_circles_tot & ...
                    avg_intensity(s-1, centers_cur, radii_cur, v) > 255*imorph.avgI_th;
                new_circles_uncertain = (new_circles_tot ...
                    & (~ new_circles_valid)) | (circle_status == -1);
                id{s-1}(new_circles_valid) = (sum(nk(1:s-2))+1):(sum(nk(1:s-2)) ...
                    + sum(new_circles_valid));
                if scan_code(1) == true
                    id{s-1}(circle_status==0) = id{s-2}(idx_cpm(circle_status==0));
                end
                id{s-1}(new_circles_uncertain) = nan;

                % extract the nuclears by using new_circles
                c{s-1} = centers_cur;
                r{s-1} = radii_cur;
                nk(s-1) = sum(new_circles_valid);
            else
                id{s-1} = [];
                warning("z resolution may be too low");
            end

            % update current focal planes
            centers_pre = centers_cur; radii_pre = radii_cur;
            centers_cur = centers_next; radii_cur = radii_next;
            s = s + 1;

            if ~isempty(bar)
                waitbar((s-1)/size(vol, 3), bar, "searching...");
            end
        end

        function sr = simr(R,r,d)
            % this function calculate the maxium similarity ratio between
            % two circle
            if d < R + r
                A = r.^2.*acos((d.^2+r.^2-R.^2)./(2*d.*r)) + ...
                    R.^2.*acos((d.^2+R.^2-r.^2)./(2*d.*R)) - ...
                    1/2*sqrt((d+r-R).*(d-r+R).*(-d+r+R).*(d+r+R));
            else
                A = 0;
            end

            sr = A./min(pi*r.^2, pi*R.^2);
        end

        function ai = avg_intensity(s, c, r, v)
            plane = double(v(:,:,s));
            ai = zeros(size(r));
            for nn = 1:length(r)
                % extract the square mean intensity
                rect = round([c(nn, :) - r(nn); c(nn, :) + r(nn)]);
                rect(1,2) = max(rect(1,2), 1);
                rect(2,2) = min(rect(2,2), size(plane, 1));
                rect(1,1) = max(rect(1,1), 1);
                rect(2,1) = min(rect(2,1), size(plane, 2));
                % correction intensity error lower than first order
                % omit the boundary rectangle weight loss
                ai(nn) = mean(plane(rect(1,2):rect(2,2),rect(1,1):rect(2,1)),"all")*4/pi ...
                    - (4/pi-1)*quantile(plane, 0.05, "all");
            end
        end

        function r = is_new(code)
            % code defination:
            % [pre_flag, next_flag, inc_flag21, inc_flag32]
            CODE_NEW = [0,1,1,-1;
                0,1,1,1;
                1,1,-1,1];
            CODE_OLD = [1,0,-1,-1;
                1,0,1,-1;
                1,1,-1,-1;
                1,1,1,-1;
                1,1,1,1];
            CODE_UNC = [0,0,1,-1];

            r = zeros(size(code, 1), 1);

            for nn = 1:length(r)
                if ismember(code(nn,:), CODE_NEW, "rows")
                    r(nn) = 1;
                elseif ismember(code(nn,:), CODE_OLD, "rows")
                    r(nn) = 0;
                elseif ismember(code(nn,:), CODE_UNC, "rows")
                    warning("low resolution or some artifact");
                    r(nn) = -1;
                else
                    error("invalid object");
                end
            end
        end
    end

    % relative-cross method
    function [v,c,r,id] = relative_cross(vol,nuclear,imopts,imorph, barflag)
        if ~exist("barflag","var")
            barflag = false;
        end

        % remap to [0, 255]
        vol = double(vol - min(vol,[],"all")) ...
            /double(max(vol,[],"all") - min(vol, [], "all"))*255;
        vol = uint8(vol);

        r_range = [round((nuclear.radius-2*nuclear.sigma)/imopts.xRes),...
            round((nuclear.radius+2*nuclear.sigma)/imopts.yRes)];

        r_range_min = min(r_range);

        if r_range_min <= 5
            % resize the volume
            vol = imresize3(vol, "Scale", [10/r_range_min, 10/r_range_min, 1]);
            imopts.xRes = imopts.xRes / (10/r_range_min);
            imopts.yRes = imopts.yRes / (10/r_range_min);
            imopts.width = size(vol, 2);
            imopts.height = size(vol, 1);
        end

        if barflag == true
            bar = waitbar(0,'preprocessing...');
        else
            bar = [];
        end
        v = vol;
        depth = size(vol, 3);

        vol = imgaussfilt3(vol);

        % use nonlinear transformation for saturation process
        % enhance the weak intensity part
        vol = double(vol)/255;
        vol = uint8(vol.^(1-imorph.gamma)*255);

        % use open operation for remove the object with radius less
        % than min radius range on each plane
        for k = 1:depth
            vol(:,:,k) = imopen(vol(:,:,k), ...
                strel('disk',r_range_min));
        end

        % using the "over-cautious and indecisive" method
        [c,r,id] = oci(v, vol, r_range, nuclear, imopts, imorph, bar);

        if barflag == true
            close(bar);
        end
    end

    function f = opfunc(mvec, vol, nuclear, imopts)
        imorph_op = struct('rehist_quantile', mvec(1), ...
                           'cir_ssy', mvec(2), ...
                           'simr_th', mvec(3), ...
                           'avgI_th', mvec(4), ...
                           'gamma',   mvec(5));
        [~,~,~,counts] = relative_cross(vol,nuclear,imopts,imorph_op);
        ivid_num = sum(cellfun(@(x)sum(isnan(x)), counts));
        subset_m = cellfun(@(x)max(x,[],"all","omitmissing"),counts, ...
            "UniformOutput",false);
        if ~all(cellfun(@(x)isempty(x), subset_m))
            vid_num = max(cell2mat(subset_m));
        else
            vid_num = 0;
        end

        % make the invalid objects less and the valid objects more
        % note that there is one-dimensional infinite square potential well
        % to punish the valid objects number
        if vid_num < search_bound(1) || vid_num > search_bound(2)
            pvid_num = 0;
        else
            potval = (vid_num-median(search_bound))^rlx_coeff;
            pvid_num = vid_num/(1+tan(pi/2*potval));
        end

        f = tradeoff_coeff*ivid_num - (1-tradeoff_coeff)*pvid_num;
    end

    % counts_view_relc: view method for relc algorithm result
    function counter_output(v,c,r,n, nuclear,imopts, varargin)
        output_ = struct('is_display', true, ...
                         'is_saving', false, ...
                         'validobj_only', false, ...
                         'format','png', ...
                         "resolution",300);
        nn = 1;
        while nn < numel(varargin)
            switch varargin{nn}
                case 'output'
                    output_ = varargin{nn+1};
            end
            nn = nn + 2;
        end
        % n is cell array, which contains idx in each cell
        if output_.is_saving
            filedir = uigetdir("select output dir");
            if isequal(filedir,0)
                output_.is_saving = false;
                disp('saving canceled');
            end
        end

        if output_.is_display || output_.is_saving
            figure;
            set(gcf, "visible", output_.is_display);
            cumn = 0;
            for k = 1:numel(n)
                imshow(v(:,:,k));
                title(gca,['slice = ',num2str(k),', #nuclear = ',...
                    num2str(sum(~isnan(n{k})))],'FontSize',15);
                if ~isempty(n{k})
                    viscircles(c{k}(~isnan(n{k}),:),r{k}(~isnan(n{k})),...
                        'LineStyle','--','Color','b','LineWidth',1);
                    if ~output_.validobj_only
                        viscircles(c{k}(isnan(n{k}),:),r{k}(isnan(n{k})),...
                            'Color','r');
                    end
                    for jj = 1:numel(r{k})
                        if ~isnan(n{k}(jj))
                            cumn = cumn + 1;
                            text(c{k}(jj,1)-0.2*nuclear.radius/imopts.xRes*(floor(log10(cumn))+1),...
                                c{k}(jj,2),num2str(n{k}(jj)),"Color",'r',...
                                "FontSize",8);
                        else
                            if ~output_.validobj_only
                                text(c{k}(jj,1)-0.5*nuclear.radius/imopts.xRes,c{k}(jj,2),...
                                    "?","Color",'y',...
                                    "FontSize",12,"FontWeight","bold");
                            end
                        end
                    end

                end
                pause(0.2);     % decay for figure refresh
                if output_.is_saving
                    exportgraphics(gca,[filedir,'\',num2str(k),'.',output_.format],...
                        "Resolution",output_.resolution);
                end
            end
            close(gcf);
        end
    end

    % rehist: re histogram the bright to enhance contrast
    function vol = rehist(vol, imopts, imorph)
        if ~exist("imorph","var")
            q = 0.97;
        else
            q =imorph.rehist_quantile;
        end 

        % sort voxel then select infsup(voxel,0.1%) as maxval
        topval = quantile(vol,0.999,"all");
        botval = quantile(vol,0.001,"all");
        vol(vol>topval) = topval;
        vol(vol<botval) = botval;

        % change the image into 8-bit image
        if imopts.dataType~="uint8"
            vol = uint8(double(vol-quantile(vol,q,"all"))*255 ...
                ./double(topval-botval));
        end
    end
end

function [lb, ub, nuclear] = genoptbound(animal_info)
marker = upper(string(animal_info.marker));
driver = upper(string(animal_info.driver));
stage = upper(string(animal_info.stage));

switch marker
    case "MCHERRY"
        switch driver
            case "OK107"
                %     req   ssy   sim   avg   gam
                lb = [0.85; 0.90; 0.70; 0.01; 0.00];
                ub = [0.99; 0.99; 1.00; 0.10; 0.75];
                switch stage
                    case "L1"
                        nuclear = struct('radius', 2.1, ...
                                         'sigma', 0.25, ...
                                         'color',"r");
                    case "L2"
                        nuclear = struct('radius', 2.1, ...
                                         'sigma', 0.25, ...
                                         'color',"r");
                    case "L3"
                        throw("Unregistered Drosophila develop stage.");
                    case "ADULT"
                        throw("Unregistered Drosophila develop stage.");
                    otherwise
                        throw("Unknown Drosophila develop stage.");
                end
            case "ORCO"
                % TODO: parameters range modified after supervised learning
                %     req   ssy   sim   avg   gam
                lb = [0.85; 0.90; 0.50; 0.02; 0.00];
                ub = [0.99; 0.99; 1.00; 0.10; 0.95];
                switch stage
                    case "L1"
                        % TODO
                    case "L2"
                        throw("Unregistered Drosophila develop stage.");
                    case "L3"
                        throw("Unregistered Drosophila develop stage.");
                    case "ADULT"
                        throw("Unregistered Drosophila develop stage.");
                    otherwise
                        throw("Unknown Drosophila develop stage.");
                end
            otherwise
                throw("Unknown driver genes on Drosophila.");
        end
    otherwise
        throw("Unknown marker genes on Drosophila.");
end
end