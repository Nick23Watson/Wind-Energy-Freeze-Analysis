%% Nicholas Watson - Energy Meteorology Final Project
% Icing on Turbines - ERA5 Reanalysis vs Observational Data
% Iowa - February 2021 
% ERA5 Reanalysis

clear; 
close all;
%% Create Modifiable Variable: WIND SPEED
windSpeedMin = 14;

%% Import Data and Read
IaERA5 = "Working\Data\Raw\IA-ERA5-FEB-1-28-2021.grib"; % Import File/Data
[A, R] = readgeoraster(IaERA5); % Read GRIB File
info = georasterinfo(IaERA5); % Extract and display the information/data
M = info.Metadata; % Save metadata table to a variable
head(M); %Display the Metadata to extract actual variable data across the month

%% Extract Data
% Temperature
idx_T = find(strcmp(M.Element,'T'));  % Extract Temperature Data Indices
% U wind
idx_U = find(strcmp(M.Element,'U')); % Extract U-component of wind Indices
% V wind
idx_V = find(strcmp(M.Element,'V')); % Extract V-Component of wind Indices
% Relative humidity
idx_RH = find(strcmp(M.Element,'R')); % Extract RH data Indices

indexTime = length(idx_T); % Set length of arrays for time indexing

% Set arrays to zeros for extracting data
[nLat, nLong, ~] = size(A(:,:,1));
T  = zeros(nLat, nLong, indexTime);
U  = zeros(nLat, nLong, indexTime);
V  = zeros(nLat, nLong, indexTime);
RH = zeros(nLat, nLong, indexTime);
timeVec = datetime(M.ValidTime(idx_T));  % time vector from metadata

%extract data and put into respective variables (T = Temp, U = U vel, V = V
%vel, RH = Rel Humid)
for k = 1:indexTime
    T(:,:,k)  = A(:,:,idx_T(k));
    U(:,:,k)  = A(:,:,idx_U(k));
    V(:,:,k)  = A(:,:,idx_V(k));
    RH(:,:,k) = A(:,:,idx_RH(k));
end

%Extract latitude and longitude bounds from the bounding box
latlim = R.LatitudeLimits;    
lonlim = R.LongitudeLimits;   

% Create vectors for each grid point
lat = linspace(latlim(1), latlim(2), R.RasterSize(1));  % rows = latitude
lon = linspace(lonlim(1), lonlim(2), R.RasterSize(2));  % cols = longitude

%% Calculations and Modifications of Data

% Calculate total wind speed (magnitude)
WindSpeed = sqrt(U.^2 + V.^2);

% Set parameters for when icing occurs: Below 0 Celsius, RH must be 90%+,
% Wind speed must be minimum - can manipulate velocity as
% different papers suggest different limits
icingFreq = (T < 2) & (RH >= 90) & (WindSpeed >= windSpeedMin);  % 3D logical array

icingHours = sum(icingFreq, 3);  % Sum over time

icingProb = icingHours / indexTime * 100;  % percentage of hours

[nLat, nLon] = size(icingHours);

% Create 2D grids for lat/lon
[lonGrid, latGrid] = meshgrid(linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), nLon), ...
                              linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), nLat));

% Flip vertically to match ERA5 orientation
icingHours = flipud(icingHours);

% Create Iowa map
figure;
ax = usamap('Iowa');  % limits automatically set for Iowa

% Apply blue-to-red colormap
colormap(jet);  
colorbar;
caxis([0 max(icingHours(:))]);  % scale colorbar to data range
title('Total Icing Hours - Iowa', 'FontSize', 14);

% Overlay Iowa border
states = shaperead('usastatelo', 'UseGeoCoords', true);  % Load US state boundaries
iowa = states(strcmp({states.Name}, 'Iowa'));            % Select Iowa
geoshow(ax, iowa, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);  % Overlay border

% % Display Final plot
% geoshow(ax, latGrid, lonGrid, icingHours, 'DisplayType', 'surface');
% colorbar;
% title(sprintf('Total Icing Hours for the Month of February - Iowa, %d m/s', windSpeedMin), 'FontSize', 14);

%% Subset for Full Month of February 
stormStart = datetime(2021,2,1,0,0,0);
stormEnd   = datetime(2021,2,28,23,59,59);

% Find indices in timeVec corresponding to storm dates
stormIdx = find(timeVec >= stormStart & timeVec <= stormEnd);

% Subset icing frequency for storm
icingFreqStorm = icingFreq(:,:,stormIdx);

% Sum icing hours over storm period
icingHoursStorm = sum(icingFreqStorm, 3);

[nLat, nLon] = size(icingHoursStorm);

% Create grids (same as before)
[lonGridStorm, latGridStorm] = meshgrid(linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), nLon), ...
                                        linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), nLat));

% Flip vertically to match ERA5 orientation
icingHoursStorm = flipud(icingHoursStorm);

%% Plot Winter Storm Icing
figure;
ax = usamap('Iowa');

% Display raster
hRasterStorm = geoshow(ax, latGridStorm, lonGridStorm, icingHoursStorm, 'DisplayType', 'surface');
alpha(hRasterStorm, 0.7);   % optional transparency
colormap(jet);
caxis([0 max(icingHoursStorm(:))]);
colorbar;
title(sprintf('Total Icing Hours (Feb 1-28), %d m/s', windSpeedMin), 'FontSize', 14);

% Add iowa borders to image
geoshow(ax, iowa, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Save figure
saveas(gcf, fullfile(pwd, sprintf('\\Working\\Results\\IA-ERA5\\IA-ERA5-%dms-FEB.png', windSpeedMin)));


