classdef NuViewer < handle
    %VIEWER This class is nuclear viewer object defination
    
    properties(Constant, Hidden)
        STATUS_SUCCESS      =   0
        STATUS_NO_HANDLE    =   -1
        REMOVED_OBJ         =   -1
        STATUS_DUPLICATE_OBJ=   -2
        
        VALID_FONT_COLOR         =   [1, 0, 0]   % valid object font color
        UNCERTAIN_FONT_COLOR     =   [1,0.5,0]   % uncertain object font color
        UNCERTAIN_OBJ_COLOR      =   [1, 1, 1]   % uncertain object color
    end

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller          % app caller, must be NuclearCounter object
        ax_caller       % axes from caller, must be uiaxes
        id_table        % NuIDTable object, for objects record
        hobj            % handle of circle objects
        cmap            % string, indicate the colormap
        volume          % data
    end

    properties(SetAccess=private, GetAccess=private)
        id              % nuclears indicator defination(identity)
        volopts         % the image options of vol
        nuclear         % nuclear setting
        center          % nuclear indicator center
        radius          % nuclear indicator radius
        color           % the colormap for viewer, m-by-3 matrix
    end

    properties(GetAccess=public, Dependent)
        obj_num_tot
        obj_max_exist
    end

    methods
        function v = get.obj_num_tot(this)
            v = this.id_table.obj_num_tot;
        end

        function v = get.obj_max_exist(this)
            v = this.id_table.obj_max_exist;
        end
    end
    
    methods(Access = public)
        function this = NuViewer(cres, vopt, animal, cmap)
            %NUVIEWER A Constructor
            arguments
                cres (1,1) struct
                vopt (1,12) table
                animal (1,1) struct
                cmap (1,1) string {mustBeMember(cmap, ...
                    ["default","hsv","jet","parula","turbo"])} = "default"
            end
            this.center = cres.center;
            this.radius = cres.radius;
            this.id = cres.id;
            this.volume = cres.volume;
            this.volopts = vopt;
            [~, ~, this.nuclear] = GenOptBound(animal);
            this.cmap = cmap;
            this.id_table = NuIDTable(this.id);
            this.gen_colormap();
        end

        function SetCaller(this, caller)
            arguments
                this (1,1) NuViewer
                caller (1,1) NuclearCounter
            end

            this.caller = caller;
            this.ax_caller = caller.UIAxes_Viewer;
        end
        
        function status = Display(this, zidx, valid_only)
            arguments
                this
                zidx (1,1) double {mustBeNonnegative, mustBeInteger} = 1
                valid_only (1,1) logical = true
            end
            %DISPLAY Disp the objects on viewer ax
            if isempty(this.ax_caller) || ~isvalid(this.ax_caller)
                status = this.STATUS_NO_HANDLE;
                return;
            end

            % update title
            this.update_title(zidx);

            id_z = this.id{zidx};
            if ~isempty(id_z)
                h_obj = cell(numel(id_z),1);
                % draw circle at the plane
                for cir_index = 1:numel(id_z)
                    cir_id = id_z(cir_index);

                    % skip the removed object
                    if cir_id == this.REMOVED_OBJ, continue; end

                    if ~isnan(cir_id)
                        self = images.roi.Circle(this.ax_caller,...
                            "Center",this.center{zidx}(cir_index,:),...
                            "Radius",this.radius{zidx}(cir_index), ...
                            "Color",this.color{zidx}(cir_index,:), ...
                            "LabelTextColor", this.VALID_FONT_COLOR, ...
                            "LineWidth",2,"MarkerSize",0.1,"Label",string(cir_id),...
                            "LabelVisible","on","FaceAlpha",0,"LabelAlpha",0,...
                            "Deletable",true,"UserData",[cir_index, zidx]);

                        % binding listener(will be auto removed if the roi was deleted)
                        addlistener(self, ...
                            'ROIMoved', @(src,~)this.roi_moved(src,[],self));

                        % binding menu: 修改标签
                        uimenu(self.ContextMenu, ...
                            "Text","修改标签",...
                            "MenuSelectedFcn",@(~,~)this.roi_modify([],[],self), ...
                            "Tag","IPTROIContextMenuModify");

                        % changed the default context menu: delete
                        self.ContextMenu.Children(end).MenuSelectedFcn ...
                            = @(~,~)this.roi_delete([],[],self);

                        h_obj{cir_index} = self;
                    else
                        if ~valid_only
                            self = images.roi.Circle(this.ax_caller,...
                            "Center",this.center{zidx}(cir_index,:),...
                            "Radius",this.radius{zidx}(cir_index), ...
                            "Color",this.UNCERTAIN_OBJ_COLOR, ...
                            "LabelTextColor", this.UNCERTAIN_FONT_COLOR, ...
                            "LineWidth",2,"MarkerSize",0.1,"Label","??",...
                            "LabelVisible","hover","FaceAlpha",0,"LabelAlpha",0,...
                            "Deletable",true, "UserData",[cir_index, zidx]);

                            % binding listener(will be auto removed if the roi was deleted)
                            addlistener(self, ...
                                'ROIMoved', @(src,~)this.roi_moved(src,[],self));

                            % binding menu: 修改标签
                            uimenu(self.ContextMenu, ...
                                "Text","修改标签",...
                                "MenuSelectedFcn",@(~,~)this.roi_modify([],[],self), ...
                                "Tag","IPTROIContextMenuModify");

                            % changed the default context menu: delete
                            self.ContextMenu.Children(end).MenuSelectedFcn ...
                                = @(~,~)this.roi_delete([],[],self);

                            h_obj{cir_index} = self;
                        end
                    end
                end
            else
                h_obj = {};
            end

            this.hobj = h_obj;
            status = this.STATUS_SUCCESS;
        end

        function status = Clear(this)
            arguments
                this
            end

            if ~isempty(this.hobj)
                for k = 1:numel(this.hobj)
                    delete(this.hobj{k}); 
                end
                this.caller.SetViewerTitle();
                status = this.STATUS_SUCCESS;
            else
                status = this.STATUS_NO_HANDLE;
            end
        end

        function status = AddObject(this, zidx, cir_id, cir_info)
            % this function add an object to data
            arguments
                this
                zidx (1,1) double {mustBePositive, mustBeInteger}
                cir_id (1,1) double % could be nan for uncertainty
                cir_info (1,3) double {mustBePositive}  % [cx,cy,r]
            end

            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end

            if ~isnan(cir_id) && ismember(cir_id, this.id{zidx})
                status = this.STATUS_DUPLICATE_OBJ;
                return;
            else
                % add new record
                this.id_table.AddObj(cir_id);
            end

            cir_index = numel(this.id{zidx});
            this.id{zidx}(cir_index+1) = cir_id;
            this.center{zidx}(cir_index+1, :) = cir_info(1:2);
            this.radius{zidx}(cir_index+1) = cir_info(3);

            this.gen_colormap();
            this.update_title(zidx);

            status = this.STATUS_SUCCESS;
        end

        function status = Refresh(this)
            % this function reorder objects identities
            % access each object and modify the identities

            % apply the refreshed identity
            [this.id, this.center, this.radius] ...
                = this.id_table.ApplyRefreshTo(this.id, ...
                                               this.center, ...
                                               this.radius, ...
                                               true);

            % regenerate the colormap
            this.gen_colormap();

            status = this.STATUS_SUCCESS;
        end

        function delete(this)
            Clear(this);

            delete(this);
        end
    end

    methods(Access = private)
        function update_title(this, zidx)
            id_z = this.id{zidx};

            % set the view title
            nuc =  sum(~isnan(id_z)) - sum(id_z==this.REMOVED_OBJ);
            unc = sum(isnan(id_z));
            title_ = sprintf("#nuclear = %d, #uncertain = %d", nuc, unc);
            this.caller.SetViewerTitle(title_);
        end

        function gen_colormap(this)
            % generate colormap, colors number as same as valid object
            % number, store as center

            % get the total objects number
            total_n = this.id_table.obj_num_tot;
            
            % generate the color map
            if this.cmap == "default"
                color_ = parula(total_n);
            else
                cfun = str2func(this.cmap);
                color_ = cfun(total_n);
            end

            % map the colors to objects
            this.color = cell(size(this.id));
            for zidx = 1:numel(this.id)
                id_z = this.id{zidx};
                for cidx = 1:numel(id_z)
                    id_z_c = id_z(cidx);
                    if ~isnan(id_z_c) && (id_z_c~=this.REMOVED_OBJ)
                        [~, id_loc] = this.id_table.IsMember(id_z_c);
                        this.color{zidx}(cidx, :) ...
                            = color_(id_loc, :);
                    end
                end
            end

        end

        function status = move_object(this, zidx, cir_index, cir_info)
            % this function modity an object in data already
            arguments
                this
                zidx (1,1) double {mustBePositive, mustBeInteger}
                cir_index (1,1) double {mustBePositive, mustBeInteger}
                cir_info (1,3) double {mustBePositive}  % [cx,cy,r]
            end
            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if cir_index > numel(this.id{zidx})
                throw(MException("NuViewer:invalidObjectIndex", ...
                    "Object index is out of range."));
            end

            % find the object
            % identity: this.id{zidx}(cir_index)
            this.center{zidx}(cir_index, :) = cir_info(1:2);
            this.radius{zidx}(cir_index, :) = cir_info(3);

            status = this.STATUS_SUCCESS;
        end

        function status = rename_object(this, zidx, cir_index, cir_id_old, cir_id_new)
            % this function rename an object in data already
            arguments
                this
                zidx (1,1) double {mustBePositive, mustBeInteger}
                cir_index (1,1) double  {mustBeNonNan} % empty for not needed: cir_id_old is not nan
                cir_id_old (1,1) double
                cir_id_new (1,1) double
            end
            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if isnan(cir_id_old) && isnan(cir_id_new)
                status = this.STATUS_SUCCESS;
                return;
            end
            if isempty(cir_index)
                if isnan(cir_id_old)
                    nanbi = isnan(this.id{zidx});
                    if sum(nanbi) > 1
                        % can't make sure the changed object
                        throw(MException("NuViewer:tooManyNaNIndex", ...
                            "Can't change label without index at this plane."));
                    else
                        % set nan -> val
                        this.id{zidx}(nanbi) = cir_id_new;
                        status = this.id_table.NaN2Val(cir_id_new);
                        this.update_title(zidx);
                        return;
                    end
                else
                    cir_index = (this.id{zidx}==cir_id_old);
                    if sum(cir_index) ~= 1
                        throw(MException("NuViewer:idNotFound", ...
                            "Identity not found."));
                    else
                        this.id{zidx}(cir_index) = cir_id_new;
                        status = this.id_table.Val2Val(cir_id_old, cir_id_new);
                        this.update_title(zidx);
                        return;
                    end
                end
            else
                if cir_index > numel(this.id{zidx})
                    throw(MException("NuViewer:invalidObjectIndex", ...
                        "Object index is out of range."));
                end
                cir_id_old = this.id{zidx}(cir_index);
                this.id{zidx}(cir_index) = cir_id_new;
                status = this.id_table.Val2Val(cir_id_old, cir_id_new);
                this.update_title(zidx);
                return;
            end
        end

        function status = remove_object(this, zidx, cir_index)
            % this function remove an object from data
            arguments
                this
                zidx (1,1) double {mustBePositive, mustBeInteger}
                cir_index (1,1) double {mustBePositive, mustBeInteger}
            end

            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if cir_index > numel(this.id{zidx})
                throw(MException("NuViewer:invalidObjectIndex", ...
                    "Object index is out of range."));
            end

            % marked the idset
            cir_id_old = this.id{zidx}(cir_index);
            status = this.id_table.DelObj(cir_id_old);

            if status == this.STATUS_SUCCESS
                this.id{zidx}(cir_index) = this.REMOVED_OBJ;
                this.update_title(zidx);
            end

        end

        function roi_moved(this, src, ~, self)
            % call Move Object
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            cinfo = [src.Center, src.Radius];

            this.move_object(zidx, cidx, cinfo);
        end

        function roi_modify(this, ~, ~, self)
            % generate a msgbox for new label acquirement
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            id_cur = this.id{zidx}(cidx);

            prompt = {'Enter new label:'};
            dlgtitle = 'rename';
            fieldsize = [1 45];
            definput = {num2str(id_cur)};

            answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
            if isempty(answer), return; end

            id_new = str2double(answer{1});

            this.rename_object(zidx, cidx, id_cur, id_new);

            % update this color
            this.gen_colormap();

            % update title
            this.update_title(zidx);

            % update object immediately
            if ~isnan(id_new)
                self.Label = string(id_new);
                self.Color = this.color{zidx}(cidx, :);
                self.LabelTextColor =  this.VALID_FONT_COLOR;
                self.LabelVisible = "on";
            else
                self.Label = "??";
                self.Color = this.UNCERTAIN_OBJ_COLOR;
                self.LabelTextColor =  this.UNCERTAIN_FONT_COLOR;
                self.LabelVisible = "hover";
            end
        end

        function roi_delete(this, ~, ~, self)
            % remove the record stored in this
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            this.remove_object(zidx, cidx);

            % remove the ROI object
            delete(self);
        end
    end

end

