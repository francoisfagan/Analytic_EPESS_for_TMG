function [phi] = sample_from_slice(slice_range) 
% Given the angle slice range, sample a point unformly from it.

    % Calculate length of intervals
    intervals = zeros(length(slice_range),1);
    for i=1:2:length(slice_range)
         intervals(i,1) = slice_range(i+1) - slice_range(i);
    end

    % Normalize
    total = sum(intervals);
    intervals = intervals/total;
    
    % Sample an interval
    index = discrete_sample(intervals,1);
    
    % Sample point in interval
    phi = rand*(slice_range(index+1) - slice_range(index)) + slice_range(index);
  
end