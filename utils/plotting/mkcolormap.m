function [colmap] = mkcolormap(maxColval,middleColval,minColval)
% [colmap] = mkcolormap(maxColval,middleColval,minColval)
% makes a one directional or two - directional colormap
% Inputs: 
%   maxColval - color for upper boundary of value range
%   middleColval - color for middle value of value range
%   minColval - color for lower boundary of value range
% Output: 
%   100x 3 (for unidirectional) or 200x 3 for diverging colormap
% (c) yulia oganian, June 2017

if nargin < 2
    middleColval = [ 1 1 1];
end

if nargin < 3
    minColval = middleColval;
end

% make sure color scale is 0-1
ColVals = [minColval; middleColval;maxColval];
for i = 1:3
    if max(ColVals(i,:)) >1
        ColVals(i,:) = ColVals(i,:)/256;
    end
end
% make colormap
if isequal(minColval, middleColval) % one-dir only
    colmap = nan(100,3);
    
    for i = 1:3
        colmap(:,i) = linspace(ColVals(2,i),ColVals(3,i), 100);
    end
else % two dir
    
    colmap = nan(200,3);
    
    for i = 1:3
        colmap(:,i) =[linspace(ColVals(1,i),ColVals(2,i), 100), linspace(ColVals(2,i),ColVals(3,i), 100)]';
    end
    
end