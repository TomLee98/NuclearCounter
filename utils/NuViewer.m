classdef NuViewer < handle
    %VIEWER This class is nuclear viewer object defination
    
    properties(SetAccess=private, GetAccess=private, Hidden)
        caller
        ax_caller
    end

    properties(SetAccess=private, GetAccess=private)
        volopts         % the image options of vol
        nuclear         % nuclear setting
        center          % nuclear indicator center
        radius          % nuclear indicator radius
        id              % nuclears indicator defination(identity)
        volume          % data
        color           % the colormap for viewer, m-by-3 matrix
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
            this.color = this.gen_colormap(cmap);
        end

        function SetCaller(this, caller)
            arguments
                this (1,1) NuViewer
                caller (1,1) NuclearCounter
            end

            this.caller = caller;
            this.ax_caller = caller.UIAxes_Viewer;
        end
        
        function [status, h_obj] = Display(this, zidx, valid_only)
            %DISPLAY Disp the objects on viewer ax
            if isempty(this.ax_caller) || ~isvalid(this.ax_caller)
                status = -1;
                return;
            end

            % set the view title
            title_ = sprintf("#nuclear = %d, #uncertain = %d", ...
                sum(~isnan(this.id{zidx})), sum(isnan(this.id{zidx})));
            this.caller.SetViewerTitle(title_);
            id_z = this.id{zidx};

            if ~isempty(id_z)
                % draw circle at the plane
                for cir_index = 1:numel(id_z)
                    if ~isnan(id_z(cir_index))
                        cir_id = this.id{zidx}(cir_index);
                        h_obj = images.roi.Circle(this.ax_caller,...
                            "Center",this.center{zidx}(cir_index,:),...
                            "Radius",rois(k, 3), ...
                            "Color",colors(rois(k, 4), :), ...
                            "LabelTextColor", round([1,1,1]-colors(rois(k, 4), :)), ... % high contrast
                            "LineWidth",1,"MarkerSize",2,"Label",labels(rois(k, 4)),...
                            "UserData",rois(k, 4),"LabelVisible","hover");



                        text(c{k}(jj,1)-0.2*nuclear.radius/imopts.xRes*(floor(log10(cir_id))+1),...
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


                


                viscircles(this.center{zidx}(~isnan(id_z),:), ...
                           this.radius{zidx}(~isnan(id_z)),...
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
    end

    methods(Access = private)
        function c = gen_colormap(this, cmap)
            % get the valid objects number
            subset_m = cellfun(@(x)max(x,[],"all","omitmissing"),this.id, ...
                "UniformOutput",false);
            if ~all(cellfun(@(x)isempty(x), subset_m))
                valid_n = max(cell2mat(subset_m));
            else
                c = [];
                return;
            end

            if cmap == "default"
                c = parula(valid_n);
            else
                cfun = str2func(cmap);
                c = cfun(valid_n);
            end
        end
    end

end

