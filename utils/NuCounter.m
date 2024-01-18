classdef NuCounter < handle
    %NUCOUNTER This class is nuclear object counter defination
    
    properties(SetAccess=private, GetAccess=private)
        % input data
        animal          % animal defination
        imorph          % input morphlogy options
        optimal         % optimal options
        nuclear         % nuclear parameters
        bkg             % camera fixed global background
    end

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller          % caller handle, caller need to implement SetProgressBar
                        % or could be [], empty for 
    end

    properties(SetAccess=private, GetAccess=private)
        % results
        center
        radius
        id
        id_table        % NuIDTable object
        volume
    end

    methods(Access=public)
        function this = NuCounter(animal, imorph, optimal, bkg)
            %COUNTER A constructor
            arguments
                animal (1,1) struct = struct('marker', "mCherry", ...
                                             'driver', "OK107", ...
                                             'stage', "L1")
                imorph (1,1) struct = struct('rehist_quantile', 0.92, ...
                                             'cir_ssy',0.95, ...
                                             'simr_th',0.85, ...
                                             'avgI_th',0.02, ...
                                             'gamma', 0.5)
                optimal (1,1) struct = struct('maxitn', 25, ...
                                              'tol', 1e-4, ...
                                              'state', 'off', ...
                                              'srange', [80,120], ...
                                              'rlx', 2, ...
                                              'trdof', 0.5)
                bkg (1,1) = 100
            end

            % generate optional boundary by animal information
            [lb, ub, nuclear] = GenOptBound(animal);
            optimal.lb = lb;
            optimal.ub = ub;

            this.animal = animal;
            this.imorph = imorph;
            this.optimal = optimal;
            this.nuclear = nuclear;
            if isnumeric(bkg) && isinf(bkg)
                this.bkg = "auto";
            else
                this.bkg = bkg;
            end
        end

        function SetCaller(this, caller)
            arguments
                this (1,1) NuCounter
                caller (1,1) NuclearCounter
            end

            this.caller = caller;
        end

        function Count(this, vol, volopts)
            % call this function for nuclears counting
            % do calculation on vol-volopts dataset
            this.caller.SetProgressBar(0);
            [center_,radius_,id_,vol_,imorph_] = count(this, vol, volopts);
            this.center = center_;
            this.radius = radius_;
            this.id = id_;
            this.volume = vol_;
            this.imorph = imorph_;
            this.id_table = NuIDTable(id_);
        end

        function r = GetResults(this)
            r = struct("center", {this.center}, ...
                       "radius", {this.radius}, ...
                       "id", {this.id}, ...
                       "volume", this.volume, ...
                       "imorph", this.imorph);
        end
    end

    methods(Access=private)
        function [center,radius,id,vol,imorph] = count(this, vol, volopts)
            % input:
            %   - this:
            %   - vol:
            %   - volopts:
            % output:
            %   - center:
            %   - radius:
            %   - id:
            %   - vol: the vol remap to [0,255]
            %   - imorph:

            % remove background at first
            if isnumeric(this.bkg)
                vol = vol - this.bkg;
            elseif isstring(this.bkg) || ischar(this.bkg)
                switch string(this.bkg)
                    case "auto"
                        % extract the angle border and calculate the median as the
                        % estimation of background
                        bkvol = vol([1:3, end-2:end], [1:3, end-2:end], :);
                        vol = vol - median(bkvol, "all");
                    otherwise
                        error("Counter:invalidBackgroundAlgorithm", ...
                            "Invalid optional argument: bkg");
                end
            end

            % rehist for enhancing contrast
            vol = this.rehist(vol, volopts);

            switch this.optimal.state
                case 'off'
                    % use the given parameters run the counter
                    [vol,center,radius,id] = ...
                        NuCounter.relative_cross(vol, volopts, this.nuclear, ...
                        this.imorph, this.caller);
                    imorph = this.imorph;
                case 'on'
                    % optimize the imaging process parameters
                    opf = @(x)NuCounter.opfunc(x, vol, volopts, this.nuclear, this.optimal);
                    % set the output function
                    gof = @(x,y,z)NuCounter.gaoutfun(x,y,z, this.caller);

                    options = optimoptions("ga",...
                        "MaxGenerations", this.optimal.maxitn, ...
                        "FunctionTolerance", this.optimal.tol, ...
                        "ConstraintTolerance", sqrt(this.optimal.tol), ...
                        "PopulationSize",35, ...    % as Drosophila cross
                        "UseParallel", true, ...
                        "MaxTime", 180, ...         % 3 min stop limit
                        "OutputFcn", gof);

                    % using genetic algorithm to solve the global optimal
                    % problem
                    if isempty(gcp("nocreate"))
                        parpool("Processes");
                    end

                    [mvx, fval] = ga(opf, numel(this.optimal.lb), [],[],[],[], ...
                        this.optimal.lb,this.optimal.ub, [], options);

                    delete(gcp("nocreate"));

                    % note that gamma as a hidden variable
                    imorph = struct('rehist_quantile', mvx(1), ...
                                    'cir_ssy',         mvx(2), ...
                                    'simr_th',         mvx(3), ...
                                    'avgI_th',         mvx(4), ...
                                    'gamma',           mvx(5));
                    this.caller.SetSystemLog(sprintf("Optimal Value: %.2f", fval));

                    % rerun the algorithm
                    [vol, center, radius, id] = ...
                        NuCounter.relative_cross(vol, volopts, this.nuclear,...
                                                 this.imorph, this.caller);
            end

            % nan_num = sum(cellfun(@(x)sum(isnan(x)), id));
            subset_m = cellfun(@(x)max(x,[],"all","omitmissing"),id, ...
                "UniformOutput",false);
            if ~all(cellfun(@(x)isempty(x), subset_m))
                valid_num = max(cell2mat(subset_m));
            else
                valid_num = 0;
            end

            this.caller.SetSystemLog(sprintf("Nuclear Objects Found: %d", valid_num));
        end

        % rehist: re histogram the bright to enhance contrast
        function vol = rehist(this, vol, volopts)
            q = this.imorph.rehist_quantile;

            % sort voxel then select infsup(voxel,0.1%) as maxval
            if numel(vol) > 1e6
                % generate 1e7 random points sample to represent the raw
                % sample
                sample_index = randi(numel(vol), 1e6, 1);
                vol_sample = vol(sample_index);
            else
                vol_sample = vol(:);
            end
            topval = quantile(vol_sample,0.999,"all");
            botval = quantile(vol_sample,0.001,"all");
            val_th = quantile(vol,q,"all");

            vol(vol>topval) = topval;
            vol(vol<botval) = botval;
            % change the image into 8-bit image
            if volopts.dataType~="uint8"
                vol = uint8(double(vol-val_th)*255./double(topval-botval));
            end
        end

        function [status, info] = anti_over_merge(this)
            % this function will anti over merge by split the except object
            % which has low Sphericity


        end

        function [status, info] = anti_over_split(this)
            % this function will anti over split by combine the except object
            % which has low Sphericity or remove object has low slices
            % number

            % try to combine the lower slices objects
        end

        function [status, info] = refine(this)
            % generate nuclears estimation and refine the results
            % 

            
        end

    end

    methods (Static)
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
                    % warning("low resolution or some artifact");
                    r(nn) = -1;
                else
                    error("invalid object");
                end
            end
        end

        % oci for relative-cross algorithm
        function [c,r,id] = oci(v, vol, r_range, nuclear, imopts, imorph, caller)
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
                            if (d_cnm(row) < radii_cur(row)) ...
                                    && (NuCounter.simr(radii_cur(row), ...
                                                      radii_next(idx_cnm(row)), ...
                                                      d_cnm(row)) >= imorph.simr_th)
                                codes(row, 2) = 1;
                            else
                                codes(row, 2) = 0;
                            end

                            % there must be codes(:,3) = 1
                            codes(row, 3) = 1;

                            % calculate the relative size
                            if (codes(row, 2) == 1) && ...
                                    (radii_next(idx_cnm(row)) >= ...
                                     radii_cur(row) - 0.5*nuclear.sigma/imopts.xRes)
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
                            if (d_cpm(row) < radii_cur(row)) ...
                                    && (NuCounter.simr(radii_cur(row), ...
                                                       radii_pre(idx_cpm(row)), ...
                                                       d_cpm(row)) >= imorph.simr_th)
                                codes(row, 1) = 1;
                            else
                                codes(row, 1) = 0;
                            end
                            % there must be codes(:,2) = 0
                            % calculate the relative size
                            if (codes(row, 1) == 0) || ...
                                    (radii_cur(row) >= ...
                                     radii_pre(idx_cpm(row)) - 0.5*nuclear.sigma/imopts.xRes)
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
                            if (d_cpm(row) < radii_cur(row)) ...
                                    && (NuCounter.simr(radii_cur(row), ...
                                                       radii_pre(idx_cpm(row)), ...
                                                       d_cpm(row)) >= imorph.simr_th) ...
                                    && (NuCounter.avg_intensity(s-2, ...
                                                                centers_pre(idx_cpm(row),:), ...
                                                                radii_pre(idx_cpm(row)), ...
                                                                v) > 255*imorph.avgI_th)
                                codes(row, 1) = 1;
                            else
                                codes(row, 1) = 0;
                            end
                            % find the next plane circle
                            if (d_cnm(row) < radii_cur(row)) ...
                                    && (NuCounter.simr(radii_cur(row), ...
                                                       radii_next(idx_cnm(row)), ...
                                                       d_cnm(row)) >= imorph.simr_th) ...
                                    && (NuCounter.avg_intensity(s, ...
                                                                centers_next(idx_cnm(row),:), ...
                                                                radii_next(idx_cnm(row)), ...
                                                                v) > 255*imorph.avgI_th)
                                codes(row, 2) = 1;
                            else
                                codes(row, 2) = 0;
                            end
                            % calculate the relative size
                            if (codes(row, 1) == 0) || ...
                                    (radii_cur(row) >= ...
                                     radii_pre(idx_cpm(row)) - 0.5*nuclear.sigma/imopts.xRes)
                                codes(row, 3) = 1;
                            else
                                codes(row, 3) = -1;
                            end
                            if (codes(row, 2) == 1) && ...
                                    (radii_next(idx_cnm(row)) >= ...
                                     radii_cur(row) - 0.5*nuclear.sigma/imopts.xRes)
                                codes(row, 4) = 1;
                            else
                                codes(row, 4) = -1;
                            end
                        end
                    end

                    % distribute the labels
                    circle_status = NuCounter.is_new(codes);
                    new_circles_tot = (circle_status==1);
                    new_circles_valid = new_circles_tot & ...
                                        (NuCounter.avg_intensity(s-1, ...
                                                                 centers_cur, ...
                                                                 radii_cur, ...
                                                                 v) > 255*imorph.avgI_th);
                    new_circles_uncertain = ...
                        (new_circles_tot & (~ new_circles_valid)) | (circle_status==-1);
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
                    if ~isempty(caller)
                        caller.SetSystemLog("Objects detected success...");
                    end
                else
                    id{s-1} = [];
                    if ~isempty(caller)
                        caller.SetSystemLog("There may be some artifacts...");
                    end
                end

                % update current focal planes
                centers_pre = centers_cur; radii_pre = radii_cur;
                centers_cur = centers_next; radii_cur = radii_next;
                s = s + 1;

                if ~isempty(caller)
                    caller.SetProgressBar((s-2)/(size(vol, 3)-2));
                end
            end
        end

        % relative-cross method
        function [v,c,r,id] = relative_cross(vol, volopts, nuclear, imorph, caller)
            if ~isempty(caller)
                caller.SetSystemLog("Preprocessing...");
            end
            % remap to [0, 255]
            vol = double(vol - min(vol,[],"all")) ...
                /double(max(vol,[],"all") - min(vol, [], "all"))*255;
            vol = uint8(vol);
            if ~isempty(caller)
                caller.SetProgressBar(0.2);
            end

            r_range = [round((nuclear.radius-2*nuclear.sigma)/volopts.xRes),...
                round((nuclear.radius+2*nuclear.sigma)/volopts.yRes)];

            r_range_min = min(r_range);

            if r_range_min <= 5
                % resize the volume
                vol = imresize3(vol, "Scale", [10/r_range_min, 10/r_range_min, 1]);
                volopts.xRes = volopts.xRes / (10/r_range_min);
                volopts.yRes = volopts.yRes / (10/r_range_min);
                volopts.width = size(vol, 2);
                volopts.height = size(vol, 1);
            end

            v = vol;
            depth = size(vol, 3);

            vol = imgaussfilt3(vol);
            if ~isempty(caller)
                caller.SetProgressBar(0.4);
            end

            % % use nonlinear transformation for saturation process
            % % enhance the weak intensity part
            vol = double(vol)/255;
            vol = uint8(vol.^(1-imorph.gamma)*255);

            % use open operation for remove the object with radius less
            % than min radius range on each plane
            for k = 1:depth
                vol(:,:,k) = imopen(vol(:,:,k), ...
                    strel('disk',r_range_min));
                if ~isempty(caller)
                    caller.SetProgressBar(0.4+k/depth*0.6);
                end
            end

            if ~isempty(caller)
                caller.SetSystemLog("Searching...");
                caller.SetProgressBar(0);
            end

            % using the "over-cautious and indecisive" method
            [c,r,id] = NuCounter.oci(v, vol, r_range, nuclear, volopts, imorph, caller);
        end

        function f = opfunc(mvec, vol, volopts, nuclear, optimal)
            imorph_op = struct('rehist_quantile', mvec(1), ...
                              'cir_ssy', mvec(2), ...
                              'simr_th', mvec(3), ...
                              'avgI_th', mvec(4), ...
                              'gamma',   mvec(5));
            [~,~,~,counts] = NuCounter.relative_cross(vol,volopts,nuclear,imorph_op, []);
            ivid_num = sum(cellfun(@(x)sum(isnan(x)), counts));
            subset_m = cellfun(@(x)max(x,[],"all","omitmissing"),counts, ...
                "UniformOutput",false);
            if ~all(cellfun(@(x)(isempty(x)||isnan(x)), subset_m))
                vid_num = max(cell2mat(subset_m));
            else
                vid_num = 0;
            end

            % make the invalid objects less and the valid objects more
            % note that there is one-dimensional infinite square potential well
            % to punish the valid objects number
            if (vid_num <  optimal.srange(1)) || (vid_num > optimal.srange(2))
                pvid_num = -inf;
            else
                potval = (vid_num-median(optimal.srange))^optimal.rlx;
                pvid_num = vid_num/(1+tan(pi/2*potval));
            end

            f = optimal.trdof*ivid_num - (1-optimal.trdof)*pvid_num;
        end

        function [state,options,optchanged] = gaoutfun(options,state,flag, caller)
            persistent parx_best      % 
            optchanged = false;

            switch flag
                case 'init'
                    % set caller info and progress bar
                    caller.SetSystemLog(sprintf("Initializing...(N=0, V=%.2f)", state.Best(end)));
                    caller.SetProgressBar(0);
                    pparx = find(state.Score==state.Best(end), 1, "last");
                    parx = state.Population(pparx, :);
                    parx_best = parx;
                    caller.SetImorphArgs(parx);
                case 'iter'
                    caller.SetSystemLog(sprintf("Evolution...(N=%d, V=%.2f)", ...
                        state.Generation, state.Best(end)));
                    caller.SetProgressBar(state.Generation/options.MaxGenerations);
                    pparx = find(state.Score==state.Best(end), 1, "last");
                    parx = state.Population(pparx, :);
                    if ~isempty(pparx)
                        parx_best = parx;
                        caller.SetImorphArgs(parx); 
                    end
                case 'done'
                    caller.SetSystemLog(sprintf("End...(N=%d, V=%.2f)", ...
                        state.Generation, state.Best(end)));
                    caller.SetProgressBar(1);
                    pparx = find(state.Score==state.Best(end), 1, "last");
                    parx = state.Population(pparx, :);
                    if ~isempty(pparx)
                        caller.SetImorphArgs(parx);
                    else
                        caller.SetImorphArgs(parx_best);
                    end
            end
        end
    end
end

