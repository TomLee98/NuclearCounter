classdef mQueue < handle
    %MQUEUE This class provide a basic queue and some operator
    
    properties (Access=private)
        length;
        data_v; % cell array is better choice
    end
    
    methods
        function this = mQueue()
            %MQUEUE Initialize an empty queue
            this.length = 0;
            this.data_v = {};
        end

        function status = push(this, d)
            this.data_v = [this.data_v, {d}];
            this.length = this.length + 1;
            status = true;
        end

        function item = pop(this)
            if ~this.isempty()
                item = this.head();
                this.data_v(1) = [];
                this.length = this.length - 1;
            else
                item = [];
            end
        end

        function item = head(this)
            if ~this.isempty()
                item = this.data_v{1};
            else
                item = [];
            end
        end

        function item = tail(this)
            if ~this.isempty()
                item = this.data_v{end};
            else
                item = [];
            end
        end

        function status = isempty(this)
            if this.size() == 0
                status = true;
            else
                status = false;
            end
        end

        function len = size(this)
            len = this.length;
        end

        function status = destroy(this)
            this.data_v = {};
            this.length = 0;
            status = true;
        end

        function q = deep_copy(this)
            q = mQueue();
            % recycle data and split the flow to new queue q
            n = this.size();
            for k = 1:n
                tmp_data = this.pop();
                q.push(tmp_data);
                this.push(tmp_data);    % append to tail
            end
        end
    end
end

