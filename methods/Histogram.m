classdef Histogram < handle
    % Histogram accumulates data in pre-defined bins
    % Histogram analyses the statistics of number of datapoints
    % in diffenent intervals (bins)

    properties(SetAccess=protected)
        nbins; % number of bins
        %what is bins?
        %the interval of Histogram. ex. [1,2);[2,3)...
        xlo; % minimum x
        xhi; % maximum x
        data; % bin contents
        underflow; % underflow counter
        overflow; % overflow counter
        %those counts datapoints outside the range (xlo,xhi)
        dx; % bin width
    end

    methods       
        % constructor
        function obj = Histogram(nbins, xlo, xhi)
            obj.nbins = nbins;
            obj.xlo = xlo;
            obj.xhi = xhi;
            obj.data = zeros(1,nbins);
            %elements of data records number of datapoints in certain bins
            %ex data[1]=30 ->30 datapoints in the first bin
            obj.underflow = 0;
            obj.overflow = 0;
            obj.dx = (xhi - xlo) / nbins;
            %overall range/number of intervals=
            %width of one interval
        end
        
        % return bin low-edges vector
        function b = bins(obj)
            bb = (obj.xlo : obj.dx : obj.xhi);
            %xlo:interval:xhi,  indeed is low-edges vector
            %ex 1,2 2,3 3,4  -> 1:1:3
            %bb is the interval vector
            %ex  bb[i] is the ith lower-edge of interval
            b = bb(1:obj.nbins);
            %why need this step? 
            %elements of b are the indexes of the interval
            %ex  b[i] is the ith interval
        end
        
        % accumulate data
        function fill(obj, v)
            %what is v?
            % the function fill v into the histogram?
            for index = 1:numel(v)
                x = v(index);
                if x < obj.xlo
                    obj.underflow = obj.underflow + 1;
                elseif x >= obj.xhi
                    obj.overflow = obj.overflow + 1;
                %counted under/overflow
                else
                    ibin = floor((x-obj.xlo)/obj.dx) + 1;
                    %what is ibin?
                    %index of the interval(bin) which contains x
                    % ex. 1,2 2,3 3,4  x=2.5, floor gives 1, ibin=2
                    obj.data(ibin) = obj.data(ibin) + 1;
                    %number of datapoints stored in the ibinth interval
                    %increments by 1
                end
            end
        end
        
        % return total number of data points accumulated
        function n = total(obj)
            n = obj.underflow + sum(obj.data) + obj.overflow;
        end
        
        % return maximum bin contents
        function n = max(obj)
            n = max(obj.data);
            %maximum number of datapoints in one bin
        end
        
        % plot the histogram
        function plot(obj)
            %plot number of datapoints against index for each bin
            figure;
            binranges = obj.bins();
            %define index vector
            bar(binranges, obj.data, 'histc');
            axis([ obj.xlo obj.xhi 0 obj.max()*1.1 ]);
            %set horizontal and vertical axis limits
        end
    end
end
