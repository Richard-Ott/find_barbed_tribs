function segment = networksegment_barbed(S,FD,DEM,A,mnratio)

%NETWORKSEGMENT_BARBED Identify river segments and compute segment geometry
%
% Syntax
%
%     segment = networksegment_barbed(S,FD,DEM,A,mnratio)
%     networksegment_barbed(S,FD,DEM,A,mnratio)
%
% Description
%
%     Object indexing the segment links between channel heads and
%     confluences, confluences and outlets and outlets and channel heads
%     
%     This function also computes the segment orientation, euclidian length
%     flow length, sinuosity, segment elevation difference, slope, drainage
%     area and strahler order
%
%     Modified to keep track of which segments share confluences (see
%     below) and change the definition of segment azimuth
%     edits from original function are highlighted
%
% Input arguments
%
%     S         stream network (STREAMobj)  
%     FD        flow direction (FLOWobj)
%     DEM       digital elevation model (GRIDobj)
%     A         flow accumulation derived from flowacc (GRIDobj)
%     mnratio   m/n ratio for steepness plots
%
% Output arguments
%
%     segment         Object indexing the segment links between channel
%                     heads and confluences, confluences and outlets and
%                     outlets and channel heads
%
%                     .IX: Global indexing refering to GRIDobj
%                     .ix: indexing refering to S (STREAMobj)
%                     .strahler: strahler order of the segment
%                     .angle: orientation of segments (0 180 degrees)
%                            flowing north -> 0° (RO edit )
%                     .length: eculidian length of segments
%                     .flength: segment length along the river network
%                     .sinuosity: segment sinuosity
%                     .dz: height difference
%                     .slope: average slope along the segment
%                     .fslope: average slope along the flow path
%                     .Amean: mean drainage area - Amean=(Amin^0.5+Amax^0.5)^2
%                     .ksn: mean steepness
%                     .confluenceInd: Index referring to segment indexes
%                     and labelling segment sharing a confluence
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,'minarea',1000);
%     segment = networksegment(S,FD,DEM,A,0.39);
%     plotsegmentgeometry(segment)
%
% 
% See also: STREAMobj, FLOWobj, streampoi, flowdistance, drainagebasins,
%           plotsegmentgeometry
% 
% Author: Philippe Steer (philippe.steer[at]univ-rennes1.fr)
% Date: 11. May, 2016

% -------------------------------------------------------------------------
% edited by Richard Ott to keep track of which networksegments share a
% common confluence.
% Also, changed the way orientation angles are computed, so that flow
% direction is considered and flowing north is set to 0°, flowing south -
% 180°
%
% 2020

%% Identify channel heads, confluences, b-confluences and outlets
Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

%% Identify basins associated to b-confluences and outlets
DB   = drainagebasins(FD,vertcat(IXbconf,IXout));
DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

%% Compute flowdistance
D = flowdistance(FD);

%% Stream order
O = streamorder(S);

%% Identify river segments
% links between channel heads and b-confluences
[~,ind11,ind12]=intersect(DBbconf,DBhead);
% links between confluences and b-confluences
[~,ind21,ind22]=intersect(DBbconf,DBconf);
% links between channel heads and outlets
[~,ind31,ind32]=intersect(DBout,DBhead);
% links between channel heads and outlets
[~,ind41,ind42]=intersect(DBout,DBconf);
% Connecting links into segments
segment.IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   segment.ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
segment.IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   segment.ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];
% Number of segments
segment.n=numel(segment.IX(:,1));

% RO edit -------
% keep track of which confluences belong together

ixc   = FD.ixc;                      % get receivers
% make fast indexing matrix
ixcix  = zeros(FD.size,'uint32');   
ixcix(FD.ix) = uint32(1):uint32(numel(FD.ix));

conf_IX  = ixc(ixcix([IXbconf(ind11)' IXbconf(ind21)']));   % get all confluence pixels, go one pixel donstream from each bcon
[~,~,ic] = unique(conf_IX);                                 % find out which confluences belong together

confluenceInds = nan(length(ibconf)/2,2); % set up output vector, save IX. Every row has the two values of the corresponding b-cons
for i = 1:length(confluenceInds)
    if length(find(ic == i)) == 2
        confluenceInds(i,:) = find(ic ==i); 
    else
        confluenceInds(i,:) = [nan, nan];  % triple/quadruple confluences cannot be used for barbed trib code, therefore skip these rare occurences
    end
end
% end EO edit ---

%% Compute segment geometry
% Compute segment strahler order
segment.strahler=double(O(segment.ix(:,1)));
% Compute segment orientation [0=E-W - 45=NE-SW 90=N-S - 135=NW-SE - 180=W-E]
x=S.x(segment.ix(:,1))-S.x(segment.ix(:,2));y=S.y(segment.ix(:,1))-S.y(segment.ix(:,2));

% RO edits ---------------------------
[theta,~] = cart2pol(S.x(segment.ix(:,1))-S.x(segment.ix(:,2)), S.y(segment.ix(:,1))-S.y(segment.ix(:,2))); % compute direction
theta = rad2deg(theta);              % convert to degreee
theta = -theta+90;                   % make sure N is up
theta(theta<0) = theta(theta<0)+360; % corrent negative values
segment.angle = theta;
% end RO edits -----------------------


% Compute segment Euclidian length
segment.length=double(sqrt((S.x(segment.ix(:,1))-S.x(segment.ix(:,2))).^2+(S.y(segment.ix(:,1))-S.y(segment.ix(:,2))).^2));
% Compute segment flow length
segment.flength=double(abs(D.Z(segment.IX(:,1))-D.Z(segment.IX(:,2))));
% Compute segment sinuosity
segment.sinuosity=double(segment.flength./segment.length);
% Segment height difference
segment.dz=double(abs(DEM.Z(segment.IX(:,1))-DEM.Z(segment.IX(:,2))));
% Segment slope along the segment
segment.slope=double(segment.dz./segment.length);
% Segment slope along the flow path
segment.fslope=double(segment.dz./segment.flength);
% Segment mean drainage area - Amean=(Amin^0.5+Amax^0.5)^2 
segment.Amean=double((sqrt(nanmin(A.Z(segment.IX(:,1)),A.Z(segment.IX(:,2))))+sqrt(nanmax(A.Z(segment.IX(:,1)),A.Z(segment.IX(:,2))))).^2);
% Segment mean steepness
segment.ksn=double(segment.fslope./(segment.Amean.*(A.cellsize^2)).^-mnratio);

% segment confluence index
segment.confluenceInd = confluenceInds;

end