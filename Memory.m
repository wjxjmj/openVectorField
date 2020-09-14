classdef Memory < handle
    properties
        data
        n
    end
    
    methods
        function self=Memory(data)
            self.data=data;
            self.n=0;
        end
    end
end