%% IDENTIFYING AND MAPPING BARBED TRIBUTARIES
% wrapper script showing the usage of the find_barb_tribs function
%
%
% K.D. Gelwick, Oct. 2022
%%
clc         % clear command window
clear       % clear workspace
close all   % close all figure windows

%% USER CHOICE -----------------------------------------------------------------------------  %
DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');  % load DEM in UTM projection
minArea = 5e5;                                  % minimum drainage area for stream initiation in m^2
mn = 0.45;                                       % set concavity (m/n ratio from Stream Power Law)
angle = 100;                                    % angle in degree above which a tributary confluence is considered barbed
min_length = 5e2;                               % minimum length of stream segment to be considered for barbed tribuataries

%% GENERATE STREAM NETWORK
disp('Generating stream network...')
FD = FLOWobj(DEM,'preprocess','carve'); % generate flow direction grid
A = flowacc(FD); % generate flow accumulation grid
S = STREAMobj(FD,'minarea',minArea,'unit','mapunits'); % generate stream network
disp('Finished generating stream network.')

%% COMPUTE RIVER SEGMENT GEOMETRY
disp('Computing network geometry...')
segment = networksegment_barbed(S,FD,DEM,A,mn); % determine stream network geometry to identify all confluences and the their angles
disp('Finished computing network geometry.')

%% IDENTIFY BARBED TRIBUTARY CONFLUENCES
disp('Identifying barbed tributaries...')
[conf_angles, barbedIX] = find_barbed_tribs(segment,angle,min_length); % identify barbed tributaries
[x,y] = ind2coord(DEM,barbedIX(:,1)); % convert linear indeces to x,y coordintes (UTM)
all_out = [barbedIX,x,y];             %  add UTM coordinates to output
barb_tribs = array2table(all_out);    % convert matrix to table
barb_tribs.Properties.VariableNames(1:4) = {'Linear_Index','Segment_Index','X_UTM','Y_UTM'}; % add variable names to table
disp('Finished identifying barbed tributaries.')

%% DETERMINE CAPTURED DRAINAGE AREA
disp('Extracting drainage area...')
drain = A.Z(barb_tribs.Linear_Index); % extract drainage area at each point based on flow accumulation grid
all_out = [all_out,drain]; % merge matrices
barb_tribs = array2table(all_out); % convert matrix to table
barb_tribs.Properties.VariableNames(1:5) = {'Linear_Index','Segment_Index','X_UTM','Y_UTM','Drainage_Area'}; % add variable names to table
disp('Finished extracting drainage area.')

%% PLOTTING
disp('Preparing map of results...')
imageschs(DEM,[],'colormap',landcolor); % display DEM with hillshade
hold on
plot(S,'k-'); % display stream network

% Plot segments
for i=1:segment.n
    plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))]...
        ,'-','Color',[.5,.5,.5],'LineWidth',2);
end

% plot barbed b confluences and scale them by drinage area of smaller
% stream entering the confluence
scat = scatter(barb_tribs.X_UTM, barb_tribs.Y_UTM,'filled','MarkerEdgeColor','k',... % display barbed tributary confluence locations
              'MarkerFaceColor','r',...
              'LineWidth',1.5);
scat.SizeData = barb_tribs.Drainage_Area/50;
% scat.SizeData = log10(barb_tribs.Drainage_Area)*10; % correlate point size to log of drainage area (linear scale is more informative in some landscapes)                

disp('Finished mapping.')