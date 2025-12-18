%% Nicholas Watson - Energy Meteorology Final Project
% Icing on Turbines - ERA5 Reanalysis vs Observational Data
% Iowa - February 2021 
% ASOS Analysis

clear; close all; clc;

%% Read CSV and Set Thresholds & Date Range
asosFile = "4187018.csv";    % set to your uploaded CSV name
outDir   = fullfile(pwd, "Working", "Results", "IA-ASOS");
if ~exist(outDir,'dir'), mkdir(outDir,'s'); end

% Icing rule chosen by user: T <= 0°C AND RH >= 90% AND
% We'll still allow AA1 codes to directly flag icing (freezing rain/drizzle)
windSpeedMin = 6;   

% Date window for analysis (edit if needed)
dateStart = datetime(2021,2,1,0,0,0);
dateEnd   = datetime(2021,2,28,23,59,59);

%% Read CSV with the data from Iowa ASOS
fprintf("Reading %s ...\n", asosFile);
T = readtable(asosFile, 'TextType','string');

% Display columns found
disp("Columns detected:");
disp(T.Properties.VariableNames);

% Required columns list
%STATION, NAME, LATITUDE, LONGITUDE, DATE, TMP, DEW, RH1, WND, AA1, CALL_SIGN
cols = T.Properties.VariableNames;

%% Parse the dates selected and filter 
if any(strcmp(cols,'DATE'))
    try
        dt_all = datetime(T.DATE, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
    catch
        dt_all = datetime(T.DATE); % fallback
    end
else
    dt_all = NaT(height(T),1);
end

% Filter rows to requested date window if dates present
maskDate = true(height(T),1);
if any(~isnat(dt_all))
    maskDate = (dt_all >= dateStart) & (dt_all <= dateEnd);
    if ~any(maskDate)
        warning('No rows fall in the requested date window; proceeding without date filtering.');
        maskDate = true(height(T),1);
    end
end
T = T(maskDate,:);
dt = dt_all(maskDate);
N = height(T);
fprintf("Rows after filtering: %d\n", N);

%% Prepare arrays for each variables
% Pre-allocate numeric vectors
tmpC = nan(N,1);
dewC = nan(N,1);
rh   = nan(N,1);
wspd = nan(N,1);
AA1present = strings(N,1);

% Helper for parsing ASOS-style subfields like "-0056,1" -> -5.6
parseASOSvalue = @(s) parseASOSSubfield_toDouble(s);

for i = 1:N
    % TMP
    if ismember('TMP', cols)
        s = T.TMP(i);
        if ~ismissing(s) && s ~=""
            v = parseASOSvalue(s);
            if ~isnan(v), tmpC(i) = v; end
        end
    end
    % DEW
    if ismember('DEW', cols)
        s = T.DEW(i);
        if ~ismissing(s) && s ~=""
            v = parseASOSvalue(s);
            if ~isnan(v), dewC(i) = v; end
        end
    end
    % RH1
    if ismember('RH1', cols)
        s = T.RH1(i);
        if ~ismissing(s) && s ~=""
            try
                p = split(strtrim(s),',');
                rv = str2double(p{1});
                if ~isnan(rv), rh(i) = rv; end
            catch
            end
        end
    end
    % WND (wind) - common ASOS format: "020,1,N,0046,1"
    if ismember('WND', cols)
        s = T.WND(i);
        if ~ismissing(s) && s ~=""
            p = split(strtrim(s),',');
            % prefer 4th element (tenths m/s) -> convert to m/s
            val = NaN;
            if numel(p) >= 4
                raw = str2double(p{4});
                if ~isnan(raw), val = raw/10; end
            end
            % fallback to 2nd element if needed
            if isnan(val) && numel(p) >= 2
                raw = str2double(p{2});
                if ~isnan(raw), val = raw/10; end
            end
            if isnan(val) && numel(p) >= 1
                raw = str2double(p{1});
                if ~isnan(raw), val = raw/10; end
            end
            wspd(i) = val;
        end
    end
    % AA1 present-weather (if available)
    if ismember('AA1', cols)
        s = T.AA1(i);
        if ismissing(s), AA1present(i) = ""; else AA1present(i) = string(s); end
    end
end

%% Compute RH if missing 
computeRH = @(Tval,Dval) (100 .* (6.112 .* exp((17.625 .* Dval) ./ (243.04 + Dval)) ./ (6.112 .* exp((17.625 .* Tval) ./ (243.04 + Tval)))));
for i = 1:N
    if isnan(rh(i)) && ~isnan(tmpC(i)) && ~isnan(dewC(i))
        rh(i) = computeRH(tmpC(i), dewC(i));
        rh(i) = min(max(rh(i),0),100);
    end
end

%% Icing flag
% Thresholds Icing if TMP <= 0°C AND RH >= 90% AND wspd > wspd_min

icing_cond = false(N,1);

hasVars = any(~isnan(tmpC)) && any(~isnan(rh)) && any(~isnan(wspd));

if hasVars
    icing_cond = (tmpC <= 0) & ...
                 (rh   >= 90) & ...
                 (wspd >= windSpeedMin);
else
    warning('TMP, RH, or wind largely missing; AA1 codes will be used where present.');
end


% Also treat certain AA1 present-weather codes as icing (override)
% Common codes: FZRA (freezing rain), FZDZ (freezing drizzle), +FZDZ, -FZRA, PL (ice pellets)
AA1_icing = false(N,1);
for i = 1:N
    s = upper(strtrim(AA1present(i)));
    if s ~= ""
        if contains(s,"FZRA") || contains(s,"FZDZ") || contains(s,"FREEZ") || contains(s,"ICE") || contains(s,"PL") || contains(s,"ICE PELLET") 
            AA1_icing(i) = true;
        end
        % also check for the AA1 '01' shorthand sometimes used
        if startsWith(s,"01")
            AA1_icing(i) = true;
        end
    end
end

icingObs = icing_cond | AA1_icing;

%% Bin hourly data
% Use station NAME for labels; fall back to CALL_SIGN then STATION
if ismember('NAME', cols)
    stationName_all = string(T.NAME);
else
    stationName_all = repmat("",N,1);
end
if ismember('CALL_SIGN', cols)
    callsign_all = string(T.CALL_SIGN);
else
    callsign_all = repmat("",N,1);
end
stationID_all = string(T.STATION);

% Lat / Lon columns from file
lat_all = nan(N,1); lon_all = nan(N,1);
if ismember('LATITUDE', cols)
    lat_all = double(T.LATITUDE);
end
if ismember('LONGITUDE', cols)
    lon_all = double(T.LONGITUDE);
end

% hourly bins
if all(~isnat(dt))
    hourBin = datetime(year(dt), month(dt), day(dt), hour(dt), 0, 0);
else
    warning('Date column not properly parsed; treating each observation separately (no hourly binning).');
    hourBin = (1:N)'; % non-hour bin fallback
end

uniqueStations = unique(stationID_all);
nStations = numel(uniqueStations);

stationLat = nan(nStations,1);
stationLon = nan(nStations,1);
stationNAME = strings(nStations,1);
stationIcingHours = zeros(nStations,1);

for s = 1:nStations
    sid = uniqueStations(s);
    mask = stationID_all == sid;
    if ~any(mask)
        continue;
    end
    % coords & name from first valid row
    idxFirst = find(mask & ~isnan(lat_all) & ~isnan(lon_all),1);
    if isempty(idxFirst)
        idxFirst = find(mask,1);
    end
    stationLat(s) = lat_all(idxFirst);
    stationLon(s) = lon_all(idxFirst);
    nm = stationName_all(idxFirst);
    if nm == "" || nm == " " || isempty(nm)
        nm = callsign_all(idxFirst);
        if nm == "" || nm==" " || isempty(nm)
            nm = sid; % fallback to station ID
        end
    end
    stationNAME(s) = nm;

    % unique hours with any icing observation for this station
    if all(~isnat(dt))
        hoursWithIcing = unique(hourBin(mask & icingObs));
        stationIcingHours(s) = numel(hoursWithIcing);
    else
        % fallback: count icing observations
        stationIcingHours(s) = sum(icingObs(mask));
    end
end

% Build results table
Results = table(uniqueStations, stationNAME, stationLat, stationLon, stationIcingHours, ...
    'VariableNames', {'StationID','NAME','LAT','LON','IcingHours'});

% Save CSV
outCSV = fullfile(outDir,'ASOS_station_icing_counts_option2.csv');
writetable(Results, outCSV);
fprintf("Wrote station icing counts to %s\n", outCSV);

%% Interpolated map generation
valid = ~isnan(Results.LAT) & ~isnan(Results.LON);
LAT = Results.LAT(valid);
LON = Results.LON(valid);
Ihours = Results.IcingHours(valid);

figure('Position',[100 100 1100 800], 'Color','w');

%% Define interpolation grid (Iowa extent)
lat_grid = 40:0.05:44.8;
lon_grid = -96.7:0.05:-90.0;
[LonG, LatG] = meshgrid(lon_grid, lat_grid);

%% Interpolate using natural neighbor
F = scatteredInterpolant(LON, LAT, Ihours, 'natural', 'none');
Igrid = F(LonG, LatG);

%% Load Iowa boundary
states = shaperead('usastatelo','UseGeoCoords',true);
IA = states(strcmp({states.Name},'Iowa'));

%% Mask outside Iowa
mask = inpolygon(LonG, LatG, IA.Lon, IA.Lat);
Igrid(~mask) = NaN;

%% Plot
ax = axesm('lambert', ...
    'MapLatLimit',[40 45], ...
    'MapLonLimit',[-97 -90], ...
    'Frame','on','Grid','on');

axis off; hold on;

% Filled interpolated field
pcolorm(LatG, LonG, Igrid);
colormap(turbo);

cb = colorbar;
ylabel(cb, 'Icing hours', 'FontSize',12);
caxis([0 max(Ihours)]);

% Overlay Iowa border
plotm(IA.Lat, IA.Lon, 'k', 'LineWidth', 2);

% Overlay ASOS station points
scatterm(LAT, LON, 60, Ihours, 'filled', 'MarkerEdgeColor','k');

% Label only highest-icing stations (top 20%)
th = prctile(Ihours,75);
label_idx = Ihours >= th;

for i = 1:height(Results)
    if ~isnan(Results.LAT(i)) && ~isnan(Results.LON(i)) && label_idx(i)
        textm(Results.LAT(i)+0.05, Results.LON(i), char(Results.NAME(i)), ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color','k');
    end
end

title({'ASOS Observed Icing Hours – 6 m/s', ...
       'February 1–28, 2021'}, ...
       'FontSize', 16);

%% Save the plot
outDirIA = fullfile(pwd, "Working", "Results", "IA-ASOS");


saveas(gcf, fullfile(outDirIA,'IA_ASOS_6ms_FEB.png'));
fprintf("Saved Iowa ASOS map\n");

%% Helper functions
function val = parseASOSSubfield_toDouble(s)
    % Parse ASOS-style values like "-0056,1" -> -5.6 (°C), or "+9999,9" -> NaN
    val = NaN;
    if ismissing(s) || s=="" 
        return
    end
    ss = strtrim(s);
    % handle JSON-like quoting or stray characters
    if startsWith(ss,'"') && endsWith(ss,'"')
        ss = ss(2:end-1);
    end
    % take first token before comma
    if contains(ss,',')
        p = split(ss,',');
        token = p{1};
    else
        token = ss;
    end
    num = str2double(token);
    if isnan(num)
        return
    end
    % missing sentinel values in ASOS are often > 5000 or = 9999
    if abs(num) > 5000
        val = NaN; return
    end
    val = num / 10;
end
