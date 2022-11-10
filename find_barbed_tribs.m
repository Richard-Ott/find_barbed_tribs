function [confluence_angles, barbedIX] = find_barbed_tribs(segment,barb_angle,min_length)
% This function maps barbed tributaries, which can be an indication for 
% river capture. 
% INPUT:    
%           - segment: structure produced by the networksegment_barbed function
%           - barb_angle: angle in degree above which a tributary
%             confluence is considered barbed
%           - min_length: minimum length of stream segment considered in
%             map units
% OUTPUT:   
%           - confluence angles (nx3): First column, value of  
%             confluence angles. Second column, IX of the b-confluence
%             with lower strahler order, 3rd column IX higher strahler
%           - barbedIX: first column, IX of all confluences with
%           barbed tribuatries, second column, index into the segment
%           structure indicating the segment with lower strahler order
%           (this can be useful for highlighting such segments in a plot)
%
% Richard Ott, 2021


n = length(segment.confluenceInd);   % number of confluences 
confluence_angles = nan(n,3);        % empty vector to store confluence angles and linear indices of confluences
j = 1;                               % counter variable
for i = 1:n
    inds = segment.confluenceInd(i,:);
    
    if isnan(inds)  % some confluence will have more than 2 tributaries and dont allow angle caclulation, skip those
        continue
    end
    
    if sum(segment.length(inds) > min_length) < 2 % check for short segments and skip those
        continue
    end
    
    % calculate confluence angle
    confluence_angles(i,1) = abs(segment.angle(inds(1)) - segment.angle(inds(2)));
    if confluence_angles(i,1) > 180   % correct when the large angle was taken instead of the smaller one 
        confluence_angles(i,1) = 360 - confluence_angles(i,1);
    end
    
    [~,lower_strahler]  = min([segment.strahler(inds(1)),segment.strahler(inds(2))]);
    [~,higher_strahler] = max([segment.strahler(inds(1)),segment.strahler(inds(2))]);
    confluence_angles(i,2) = segment.IX(inds(lower_strahler));   % store first lower strahler bconfluence coordinate
    confluence_angles(i,3) = segment.IX(inds(higher_strahler));  % store higher strahler bconfluence coordinate
    
    barbed = confluence_angles(i,1) > barb_angle;    % check if the confluence angle is larger than the one indicated
    if barbed
        barbedIX(j,1) = segment.IX(inds(lower_strahler));  % if the tributary is barbed, store IX
        barbedIX(j,2) = inds(lower_strahler);              % store the segment index of the lower strahler order stream for plotting 
        j = j+1; % add 1 to counter
    end
end


end

