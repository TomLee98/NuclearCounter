classdef mStack < handle
    %MSTACK This class provide a basic stack and some operator
    
    properties (Access=private)
        length;
        data_v; % cell array is better choice
    end
    
    methods
        function this = mStack()
            %MSTACK Initialize an empty stack
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
                item = this.top();
                this.data_v(end) = [];
                this.length = this.length - 1;
            else
                item = [];
            end
        end

        function item = top(this)
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

        function s = deep_copy(this)
            s = mStack();
            t = mStack();
            n = this.size();
            for k = 1:n
                tmp_data = this.pop();
                s.push(tmp_data);
                t.push(tmp_data);
            end
            s.inverse();
            % recover this
            while ~t.isempty()
                this.push(t.pop());
            end
        end
    end

    methods(Access = private)
        % this function help for deep copy
        % which can inverse data sequence in place abstractly
        function this = inverse(this)
            this.data_v = fliplr(this.data_v);
        end
    end
end