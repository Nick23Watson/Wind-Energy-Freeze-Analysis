%% Nicholas Watson - Energy Meteorology Final Project
% ERA5 Single-Point Icing Analysis for Specific Locations
% Iowa - February 2021

clear; 
close all; 
clc;

%% Settings

% ERA5 file
TxERA5 = "Working\Data\Raw\IA-ERA5-FEB-1-28-2021.grib";

% Station location (match ASOS station)
stationName = "Ottumwa Industrial Airport (KOTM)";
stationLat  = 41.10079;
stationLon  = -92.4445;

% Date range (Edit for desired times)
dateStart = datetime(2021,2,1,0,0,0);
dateEnd   = datetime(2021,2,28,23,59,59);

% Icing thresholds (MATCH ASOS)
T_thresh  = 0;     % °C
RH_thresh = 90;    % %
windMin   = 2;     % m/s

%% read grib file

[A, R] = readgeoraster(TxERA5);
info = georasterinfo(TxERA5);
M = info.Metadata;

%% extract the relevant variables from the grib file

idx_T  = find(strcmp(M.Element,'T'));
idx_U  = find(strcmp(M.Element,'U'));
idx_V  = find(strcmp(M.Element,'V'));
idx_RH = find(strcmp(M.Element,'R'));

timeVecAll = datetime(M.ValidTime(idx_T));

[nLat, nLon, ~] = size(A(:,:,1));
nTime = numel(idx_T);

T  = zeros(nLat,nLon,nTime);
U  = zeros(nLat,nLon,nTime);
V  = zeros(nLat,nLon,nTime);
RH = zeros(nLat,nLon,nTime);

for k = 1:nTime
    T(:,:,k)  = A(:,:,idx_T(k));
    U(:,:,k)  = A(:,:,idx_U(k));
    V(:,:,k)  = A(:,:,idx_V(k));
    RH(:,:,k) = A(:,:,idx_RH(k));
end

%% Find the nearest grid point to the coordinates of the selected station

lat = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1));
lon = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2));

[~, iLat] = min(abs(lat - stationLat));
[~, iLon] = min(abs(lon - stationLon));

fprintf("Using ERA5 grid point: %.2f°N, %.2f°E\n", lat(iLat), lon(iLon));

%% Set variables in arrays

T_pt  = squeeze(T(iLat,iLon,:));
U_pt  = squeeze(U(iLat,iLon,:));
V_pt  = squeeze(V(iLat,iLon,:));
RH_pt = squeeze(RH(iLat,iLon,:));

WSP_pt = sqrt(U_pt.^2 + V_pt.^2);

%% Set time for data to be observed from for the plots

maskTime = (timeVecAll >= dateStart) & (timeVecAll <= dateEnd);

dt  = timeVecAll(maskTime);
TMP = T_pt(maskTime);
RH  = RH_pt(maskTime);
WSP = WSP_pt(maskTime);

RH(RH > 100) = 100;   % ERA5 sometimes gets RH too high, correct


fprintf("ERA5 points used: %d\n", numel(dt));

%% Logic for icing using thresholds

icing = (TMP <= T_thresh) & ...
        (RH  >= RH_thresh) & ...
        (WSP >= windMin);

%% Place hourly data into bins to be references in plots

hourBin     = datetime(year(dt),month(dt),day(dt),hour(dt),0,0);
uniqueHours = unique(hourBin);

hourlyIcing = zeros(numel(uniqueHours),1);
icingHours  = 0;

for h = 1:numel(uniqueHours)
    idx = hourBin == uniqueHours(h);
    if any(icing(idx))
        hourlyIcing(h) = 1;
        icingHours = icingHours + 1;
    end
end

fprintf("Total ERA5 icing hours: %d\n", icingHours);

%% Plots

figure('Position',[200 100 1100 900],'Color','w');
tiledlayout(4,1,'TileSpacing','compact');

% Temperature Plot
nexttile
plot(dt,TMP,'k','LineWidth',1.3); hold on;
yline(T_thresh,'r--','HandleVisibility','off');
text(max(dt), T_thresh,'  0°C', ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
ylabel('Temp (°C)');
title(sprintf('%s (ERA5 Grid Point)', stationName));
grid on;

% Relative Humidty plot
nexttile
plot(dt,RH,'b','LineWidth',1.3); hold on;
yline(RH_thresh,'r--','HandleVisibility','off');
text(max(dt), RH_thresh,'  90%', ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
ylabel('RH (%)');
grid on;

% Plot wind speed
nexttile
plot(dt,WSP,'Color',[0.2 0.6 0.2],'LineWidth',1.3); hold on;
yline(windMin,'r--','HandleVisibility','off');
text(max(dt), windMin, sprintf('  %g m/s',windMin), ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
ylabel('Wind (m/s)');
grid on;

% Icing hours plot
nexttile
stairs(uniqueHours,hourlyIcing,'LineWidth',1.6);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'No Icing','Icing'});
ylabel('Icing');
xlabel('Time');
grid on;

sgtitle(sprintf(['ERA5 Single-Point Icing Analysis\n', ...
    '%s | %.2f°N %.2f°E\n', ...
    '%s to %s | Total Icing Hours = %d'], ...
    stationName, stationLat, stationLon, ...
    datestr(dateStart), datestr(dateEnd), icingHours), ...
    'FontWeight','bold');

%% save the plot

outDir = fullfile(pwd,'Working','Results','ERA5-Regional-IA');
if ~exist(outDir,'dir'), mkdir(outDir); end

saveas(gcf, fullfile(outDir, ...
    sprintf('ERA5_%s.png', replace(stationName,' ','_'))));

fprintf("ERA5 plot saved successfully.\n");
