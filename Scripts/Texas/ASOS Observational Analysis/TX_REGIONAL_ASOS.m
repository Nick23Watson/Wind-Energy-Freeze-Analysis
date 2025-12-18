%% Nicholas Watson - Energy Meteorology Final Project
% ASOS Observations – Single Station Analysis

clear; 
close all; 
clc;

%% Set thresholds and dates

asosFile  = "4186137.csv";        % ASOS CSV
stationID = 72255012912;          % NUMERIC station ID (e.g., KMAF)

% Date range (EDIT)
dateStart = datetime(2021,2,14,0,0,0);
dateEnd   = datetime(2021,2,17,23,59,59);

% Icing thresholds
T_thresh  = 0;        % °C
RH_thresh = 90;       % %
windMin   = 2;        % m/s

%% Read the data

fprintf("Reading ASOS file...\n");
T = readtable(asosFile,'TextType','string');

% Filter ID of weather stations
if isnumeric(T.STATION)
    maskStation = T.STATION == stationID;
else
    maskStation = string(T.STATION) == string(stationID);
end

T = T(maskStation,:);

if isempty(T)
    error("No data found for station ID %s", string(stationID));
end

stationName = unique(T.NAME);
stationLat  = T.LATITUDE(1);
stationLon  = T.LONGITUDE(1);

fprintf("Station: %s\n", stationName);

%% Select dates and times

dt = datetime(T.DATE, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss');

maskTime = (dt >= dateStart) & (dt <= dateEnd);
T  = T(maskTime,:);
dt = dt(maskTime);

fprintf("Observations used: %d\n", height(T));

%% Parse through each of the variables and save the data 

N = height(T);

TMP = nan(N,1);
DEW = nan(N,1);
RH  = nan(N,1);
WSP = nan(N,1);

for i = 1:N

    % Temperature
    if ~ismissing(T.TMP(i)) && T.TMP(i) ~= ""
        TMP(i) = parseASOS(T.TMP(i));
    end

    % Dewpoint
    if ~ismissing(T.DEW(i)) && T.DEW(i) ~= ""
        DEW(i) = parseASOS(T.DEW(i));
    end

    % -RH
    if ~ismissing(T.RH1(i)) && T.RH1(i) ~= ""
        p = split(T.RH1(i),',');
        RH(i) = str2double(p{1});
    end

    % -Wind Speed
    if ~ismissing(T.WND(i)) && T.WND(i) ~= ""
        p = split(T.WND(i),',');
    
        if numel(p) >= 4
            raw = str2double(p{4});
    
            % Reject missing/sentinel values
            if ~isnan(raw) && raw > 0 && raw < 600 % make sure no absurd wind speeds are used
                % raw is tenths of m/s
                WSP(i) = raw / 10;
            else
                WSP(i) = NaN;
            end
        end
    end

end

%% Need to calculate RH from T/Dew

for i = 1:N
    if isnan(RH(i)) && ~isnan(TMP(i)) && ~isnan(DEW(i))
        RH(i) = 100 * ...
            exp((17.625*DEW(i))/(243.04+DEW(i))) / ...
            exp((17.625*TMP(i))/(243.04+TMP(i)));
    end
end

%% Save data for icing hours based on observed data and parameters

icing = (TMP <= T_thresh) & (RH >= RH_thresh);

if any(~isnan(WSP))
    icing = icing & (WSP >= windMin);
end

%% Save each of the data from each our in bins
hourBin = datetime(year(dt),month(dt),day(dt),hour(dt),0,0);
uniqueHours = unique(hourBin);

icingHours = 0;
hourlyIcing = zeros(numel(uniqueHours),1);

for h = 1:numel(uniqueHours)
    idx = hourBin == uniqueHours(h);
    if any(icing(idx))
        icingHours = icingHours + 1;
        hourlyIcing(h) = 1;
    end
end

fprintf("Total icing hours: %d\n", icingHours);

%% Generate plot with each of the variables across the time range

figure('Position',[200 100 1100 900],'Color','w');
tiledlayout(4,1,'TileSpacing','compact');

% Make Temperature section of Plot
nexttile
plot(dt,TMP,'k','LineWidth',1.3); hold on;
yline(T_thresh,'r--','HandleVisibility','off');
text(max(dt), T_thresh, '  0°C', ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','left');
ylabel('Temp (°C)');
title(sprintf('%s (Station %d)', stationName, stationID));

% Make Relative Humidty Section of Plot
nexttile
plot(dt,RH,'b','LineWidth',1.3); hold on;
yline(RH_thresh,'r--','HandleVisibility','off');
text(max(dt), RH_thresh, '  90%', ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','left');
grid on;
ylabel('RH (%)');

% Make Wind Speed Plot
nexttile
plot(dt,WSP,'Color',[0.2 0.6 0.2],'LineWidth',1.3); hold on;
yline(windMin,'r--','HandleVisibility','off');
text(max(dt), windMin, sprintf('  %g m/s',windMin), ...
    'Color','r','FontWeight','bold', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','left');
grid on;
ylabel('Wind (m/s)');

% Insert section for icing at bottom of plot
nexttile
stairs(uniqueHours,hourlyIcing,'LineWidth',1.6);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'No Icing','Icing'});
grid on;
xlabel('Time');
ylabel('Icing');

sgtitle(sprintf(['ASOS Single-Station Icing Analysis\n', ...
    '%s | %.2f°N %.2f°E\n', ...
    '%s to %s | Total Icing Hours = %d'], ...
    stationName, stationLat, stationLon, ...
    datestr(dateStart), datestr(dateEnd), icingHours), ...
    'FontWeight','bold');

%% Save the image 

outDir = fullfile(pwd,'Working','Results','ASOS-Regional-TX');
if ~exist(outDir,'dir'), mkdir(outDir); end

saveas(gcf, fullfile(outDir, ...
    sprintf('ASOS_%d.png', stationID)));

fprintf("Saved plot successfully.\n");

%% Include Helper Files

function val = parseASOS(s)
    val = NaN;
    if ismissing(s) || s==""
        return
    end
    p = split(s,',');
    num = str2double(p{1});
    if abs(num) > 5000
        val = NaN;
    else
        val = num/10;
    end
end
