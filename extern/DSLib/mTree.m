classdef mTree < handle
    %MTREE This class provide a basic tree and some operator
    
    properties(Access = private)
        tree_v      % cell array for mTree Node
        id_table    % array of sub nodes id
        data        % any type of data
        length      % size of children
    end

    properties(GetAccess=public, SetAccess=private)
        id          % the unique identity of a node in one tree
    end

    properties(Dependent, GetAccess=public)
        isleaf      %
    end
    
    methods
        % set the identity only once
        function setupImpl(this, id)
            this.id = id;
        end

        function this = mTree(id, data)
            arguments
                id (1,1) double {mustBeNonnegative, mustBeInteger} = 0
                data  = []
            end 
            %MBTREE Initialize an empty tree
            this.length = 0;
            this.data = data;
            this.tree_v = {};
            this.id_table = [];
            % Init identity only once
            setupImpl(this, id);
        end
        
        function status = addnode(this, node)
            if this.is_node_exist(node.id)
                % make sure any node of tree d
                this.tree_v = [this.tree_v, {node}];
                this.length = this.length + 1;

                status = true;
            else
                status = false;
            end
        end

        function status = rmnode_by_val(this, val)
            % remove node by its value, maybe some nodes have the same data
            marker = false(1, this.length);
            for k = 1:this.length
                if this.tree_v{k}.getdata() == val
                    marker(k) = true;
                end
            end
            
            if any(marker)
                this.tree_v(marker) = [];
                this.length = this.length - sum(marker);

                status = true;
            else
                status = false;
            end
        end

        function status = rmnode_by_id(this, id)
            % remove node by its id, which is unique in a tree
            for k = 1:this.length
                if this.tree_v{k}.id == id
                   this.tree_v(k) = [];
                   this.length = this.length - 1;

                   status = true;
                   return;
                end
            end

            status = false;
        end

        function node = getnode_by_id(this, id)
            % get node by its id
            for k = 1:this.length
                if this.tree_v{k}.id == id
                   node = this.tree_v{k};
                   return;
                end
            end
            node = [];
        end

        function nodes = getnode_by_val(this, val)
            % return a n-by-1 mTree array, each element has data equal to
            % val
            nodes = [];
            for k = 1:this.length
                if this.tree_v{k}.getdata() == val
                   nodes = [nodes, this.tree_v{k}]; %#ok<AGROW>
                end
            end
        end

        function sz = numel(this)
            sz = this.length;
        end

        function status = setdata(this, d)
            this.data = d;

            status = true;
        end

        function d = getdata(this)
            d = this.data;
        end

        function tf = get.isleaf(this)
            if numel(this) == 0
                tf = true;
            else
                tf = false;
            end
        end

        function tf = eq(lhs, rhs)
            % this function compare lhs and rhs are same tree by compare
            % the identity
            tf = (lhs.id == rhs.id);
        end

        function tf = veq(lhs, rhs)
            % this function compare lhs and rhs are trees which contains 
            % the same data
            tf = (lhs.getdata() == rhs.getdata());
        end

        function tree = copy(this)
            % this function is deep copy of this tree
            % which has same id and data
            
        end
    end

    methods (Access = private)
        function tf = is_node_exist(this, id)
            % find if the tree has node id equal to the same id
            if (this.id ~= id) && ~ismember(id, this.id_table)
                tf = false;
            else
                tf = true;
            end
        end
    end

    methods(Static)
        
    end
end

