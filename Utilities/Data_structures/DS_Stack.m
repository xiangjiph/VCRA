classdef DS_Stack < handle
    properties(Access = private)
       data_cell = {};
       pointer_idx = 0;
    end
%% Static method    
    methods(Static)
        
    end
%% Non-static method
    methods
%         function obj = DS_Stack(obj)
%             obj.data_cell = {};
%             obj.pointer_idx = 0;
%         end        
        function obj = push(obj, new_data)
            obj.pointer_idx = obj.pointer_idx + 1;
            obj.data_cell{obj.pointer_idx} = new_data;
        end
        
        function output = pop(obj)
            if obj.pointer_idx > 0
                output = obj.data_cell{obj.pointer_idx};
                obj.data_cell(obj.pointer_idx) = [];
                obj.pointer_idx = obj.pointer_idx - 1;
            else
                error('Empty stack');
            end
        end
        
        function output = isEmpty(obj)
            output = (obj.pointer_idx == 0);
        end
        
        function output = numElem(obj)
            output = obj.pointer_idx;
        end        
        
        function obj = clearAll(obj)
            obj.data_cell = {};
            obj.pointer_idx = 0;
        end
    end 

end