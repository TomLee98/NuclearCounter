classdef NuIDTable < handle
    %NUIDTABLE This class is helper for NuViewer/NuCounter as idtable
    % defination
    
    properties(SetAccess = private, GetAccess = private)
        id_table        % n-by-2 array for data recording
        nan_num         % 1-by-1 nonnegtive integer, nan recording
        rm_num          % 1-by-1 nonnegtive integer, removed recording
    end

    properties(GetAccess=public, Dependent)
        obj_num_tot     % total objects number(include objs with removed tag)
        obj_num_exist   % only objects without removed tag number
        obj_max_exist   % maximum identity of exist objects
    end

    properties(Constant, Hidden)
        STATUS_SUCCESS      =   0
        STATUS_NO_OBJ       =   -1
        STATUS_INF_OBJ      =   -2
        STATUS_NAN_OBJ      =   -3
        STATUS_NO_REFRESH   =   -4

        REMOVED_OBJ = -1
    end

    methods
        function v = get.obj_num_tot(this)
            v = size(this.id_table, 1);
        end

        function v = get.obj_num_exist(this)
            exist_loc = (this.id_table(:,2)~=0);
            v = sum(exist_loc, "all");
        end

        function v = get.obj_max_exist(this)
            exist_loc = (this.id_table(:,2)~=0);
            v = max(this.id_table(exist_loc, 1),[],"all");
        end
    end    

    methods(Access=public)
        function this = NuIDTable(idset)
            %NUIDTABLE A constructor
            arguments
                idset (:,1) cell
            end

            idtable_ = cell2mat(idset);
            this.nan_num = sum(isnan(idtable_),"all");
            this.rm_num = 0;

            idtable_c1 = unique(idtable_);
            idtable_c1(isnan(idtable_c1)) = [];

            idtable_c2 = histcounts(cell2mat(idset),"BinWidth",1, ...
                "BinLimits",[1, max(idtable_c1,[],"all","omitmissing")+1]);

            this.id_table = [idtable_c1, idtable_c2'];
        end
        
        function status = DelObj(this, id)
            if isnan(id)
                if this.nan_num > 0
                    this.nan_num = this.nan_num - 1;
                    this.rm_num = this.rm_num + 1;
                    status = this.STATUS_SUCCESS;
                    return; 
                else
                    status = this.STATUS_NO_OBJ;
                    return;
                end
            else
                deloc = (this.id_table(:,1)==id);
                if any(deloc)
                    if this.id_table(deloc, 2) > 0
                        this.id_table(deloc, 2) ...
                            = this.id_table(deloc, 2) - 1;
                        this.rm_num = this.rm_num + 1;
                        status = this.STATUS_SUCCESS;
                        return;
                    else
                        status = this.STATUS_NO_OBJ;
                        return;
                    end
                else
                    status = this.STATUS_NO_OBJ;
                    return;
                end
            end
        end

        function status = AddObj(this, id)
            if isnan(id)
                this.nan_num = this.nan_num + 1;
            else
                if isinf(id)
                    status = this.STATUS_INF_OBJ;
                    return;
                end

                adloc = (this.id_table(:,1)==id);
                if any(adloc)
                    this.id_table(adloc, 2) ...
                        = this.id_table(adloc, 2) + 1;
                else
                    % add new item
                    this.id_table = [this.id_table; [id, 1]];
                end
            end

            status = this.STATUS_SUCCESS;
        end

        function status = Val2Val(this, id1, id2)
            % change from id1 to id2
            % <=> delete id1 and add id2
            status = this.DelObj(id1);

            if status == this.STATUS_SUCCESS
                status = this.AddObj(id2);
            end
        end

        function status = NaN2Val(this, id)
            status = this.Val2Val(nan, id);
        end

        function status = Val2NaN(this, id)
            status = this.Val2Val(id, nan);
        end

        function [tf, pos] = IsMember(this, id)
            adloc = (this.id_table(:,1)==id);
            tf = any(adloc);
            if nargout == 2
                pos = find(adloc);
            end
        end

        function [id, c, r] = ApplyRefreshTo(this, id, c, r, reorder)
            [status, idmap] = this.refresh(reorder);

            if status == this.STATUS_NO_REFRESH, return; end

            for zidx = 1:numel(id)
                % skip empty
                if isempty(id{zidx}), continue; end

                % remove objects
                obj_loc_del = (id{zidx}==this.REMOVED_OBJ);
                if any(obj_loc_del)
                    id{zidx}(obj_loc_del) = [];
                    c{zidx}(obj_loc_del, :) = [];
                    r{zidx}(obj_loc_del) = [];
                end

                % remap label
                if reorder == true
                    for cidx = 1:numel(id{zidx})
                        if ~isnan(id{zidx}(cidx))
                            cir_loc = (idmap(:,1)==id{zidx}(cidx));
                            id{zidx}(cidx) = idmap(cir_loc, 2);
                        end
                    end
                end

            end
        end
    end

    methods(Access=private)
        function [status, remap] = refresh(this, reorder)
            arguments
                this
                reorder (1,1) logical = false
            end
            % refresh will 
            % 1. remove objects with zero record
            % 2. reorder the record if true
            % 3. reset rm_num flag to zero

            rm_loc = (this.id_table(:,2)==0);
            if ~any(rm_loc)
                status = this.STATUS_NO_REFRESH;
                remap = [];
                return;
            else
                this.id_table(rm_loc, :) = [];
                if reorder == true
                    remap = [this.id_table(:,1), (1:size(this.id_table, 1))'];
                    this.id_table(:,1) = (1:size(this.id_table, 1))';
                else
                    remap = repmat(this.id_table(:,1),1,2);
                end
            end

            this.rm_num = 0;

            status = this.STATUS_SUCCESS;
        end
    end
end