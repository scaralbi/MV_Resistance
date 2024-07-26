%% Import Data
cd ../Data/Spectroscopy

All = readmatrix('20072024_15hour_assay.xlsx', 'Sheet', 'All');


cd ../../Scripts




%% Import Echem


cd ../Data/Echem_MV

DPV = readmatrix('AS_SWV_DPV_0hr15hr.xlsx', 'Sheet', 'DPV');
SWV = readmatrix('AS_SWV_DPV_0hr15hr.xlsx', 'Sheet', 'SWV');




cd ../../Scripts

%% Params
cmap = brewermap(6,'Blues');
cmap_blanked = cmap(2:end, :);
cmapT = brewermap(3,'Set2');

% Seaborn Tab10 palette

% Seaborn Tab10 palette
cmapw(1,:) = [72 118 177]./255;
cmaps(2,:) = [227 126 35]./255;
cmaps(3,:) = [97 159 58]./255;
cmaps(4,:) = [184 43 44]./255;
cmaps(5,:) = [136 103 185]./255;
cmapw(6,:) = [127 87 77]./255;
cmaps(7,:) = [202 118 190]./255;
cmaps(8,:) = [126 126 126]./255;
cmaps(9,:) = [187 188 56]./255;
cmaps(10,:) = [109 189 205]./255;


DPV_V = DPV(:,1)+0.2; %convert to SHE;
DPV_A = DPV(:,2:end);
DPV_strains = DPV_A(:, 2:end);
DPV_blank = repmat(DPV_A(:,1), 1, size(DPV_strains,2));

SWV_V = SWV(:,1)+0.2; %convert to SHE
SWV_A = SWV(:,2:end);
SWV_strains = SWV_A(:, 2:end);
SWV_blank = repmat(SWV_A(:,1), 1, size(SWV_strains,2));

%blank measurements
SWV_blanked = SWV_strains - SWV_blank;
DPV_blanked = DPV_strains - DPV_blank;

k_dp = find((-0.5<=DPV_V) & (DPV_V<=-0.2));
k_sw = find((-0.5<=SWV_V) & (SWV_V<=-0.2));

k_sw0 = find((-0.6<=SWV_V) & (SWV_V<=0.1));
k_sw15 = find((-0.67<=SWV_V) & (SWV_V<=0.03));


DPV_area_V = DPV_V(k_dp(1):k_dp(end), 1);
SWV_area_V = SWV_V(k_sw(1):k_sw(end), 1);

SWV_area_V0 = SWV_V(k_sw0(1):k_sw0(end), 1);
SWV_area_V15 = SWV_V(k_sw15(1):k_sw15(end), 1);

SWV0 = [SWV_blanked(:,1), SWV_blanked(:,3), SWV_blanked(:,5), SWV_blanked(:,7), SWV_blanked(:,9), SWV_blanked(:,11), SWV_blanked(:,13), SWV_blanked(:,15), SWV_blanked(:,17), SWV_blanked(:,19)];
SWV15 = [SWV_blanked(:,2), SWV_blanked(:,4), SWV_blanked(:,6), SWV_blanked(:,8), SWV_blanked(:,10), SWV_blanked(:,12), SWV_blanked(:,14), SWV_blanked(:,16), SWV_blanked(:,18), SWV_blanked(:,20)];

DPV0 = [DPV_blanked(:,1), DPV_blanked(:,3), DPV_blanked(:,5), DPV_blanked(:,7), DPV_blanked(:,9), DPV_blanked(:,11), DPV_blanked(:,13), DPV_blanked(:,15), DPV_blanked(:,17), DPV_blanked(:,19)];
DPV15 = [DPV_blanked(:,2), DPV_blanked(:,4), DPV_blanked(:,6), DPV_blanked(:,8), DPV_blanked(:,10), DPV_blanked(:,12), DPV_blanked(:,14), DPV_blanked(:,16), DPV_blanked(:,18), DPV_blanked(:,20)];


DPV_area = DPV_blanked(k_dp(1):k_dp(end), :);
SWV_area = SWV_blanked(k_sw(1):k_sw(end), :);

SWV0_area = SWV0(k_sw0(1):k_sw0(end), :);
SWV15_area = SWV15(k_sw15(1):k_sw15(end), :);


DPV0_area = DPV0(k_dp(1):k_dp(end), :);
DPV15_area = DPV15(k_dp(1):k_dp(end), :);

DPV_int = trapz(max(DPV_area,0));

SWV_int = trapz(max(SWV_area,0));

SWV0_int = trapz(max(SWV0_area,0));
SWV15_int = trapz(max(SWV15_area,0));

DPV0_int = trapz(max(DPV0_area,0));
DPV15_int = trapz(max(DPV15_area,0));




DPV_Area = [DPV0_int; DPV15_int];
SWV_Area = [SWV0_int; SWV15_int];


%SWV calibration from peak area from calibration experiments
% y(area) = 3.40x(mv conc) + 4.82
% x = (y - 4.82)/3.40


%SWV calibration from peak area
% y(area) = 7.05x(mv conc) + 2.13
% x = (y - 2.13) / 7.05



SW_MV0 = (SWV0_int-4.82)./3.40;
SW_MV15 = (SWV15_int-4.82)./3.40;
SWMV = [SW_MV0; SW_MV15];
SW_diff = (SW_MV15./SW_MV0).*100; %  fraction of the initial MV that remains in the supernatant after 15 hours expressed as percentage


DP_MV0 = (DPV0_int-2.13)./7.05;
DP_MV15 = (DPV15_int-2.13)./7.05;
DPMV = [DP_MV0; DP_MV15];
DP_diff = (DP_MV15./DP_MV0).*100; %  fraction of the initial MV that remains in the supernatant after 15 hours expressed as percentage


Diff = [DP_diff; SW_diff];


%%  Assay


titles = {'wt_Nixon','mvR01_Nixon', 'mvR02_Nixon', 'mvR03_Nixon', 'mvR06_Nixon', 'wt_Howe', 'mvR09_Howe', 'mvR10_Howe', 'mvR11_Howe', 'mvR12_Howe'};


OD750 = All(1:3, 3:end);

OD680 = All(4:6, 3:end);

DCF = All(10:12, 3:end);


DCF_norm = DCF./OD750;

DCF_avg = mean(DCF_norm,1);
DCF_std = std(DCF_norm,0,1);




%% Function to find peaks around a specified x value
function [peak_position, peak_height] = find_peak_around_target(voltage, current, target_x, tolerance)
    % Find peaks
    [peak_heights, peak_indices] = findpeaks(current, 'MinPeakProminence', 0.2); % Adjust prominence threshold as needed
    peak_positions = voltage(peak_indices);
    
    % Filter peaks around the target_x within the tolerance range
    within_tolerance = abs(peak_positions - target_x) <= tolerance;
    peak_positions = peak_positions(within_tolerance);
    peak_heights = peak_heights(within_tolerance);
    
    % If no peaks found, return y value at target_x
    if isempty(peak_positions)
        [~, closest_idx] = min(abs(voltage - target_x));
        peak_position = voltage(closest_idx);
        peak_height = current(closest_idx);
    else
        % Find the highest peak within the range
        [~, max_idx] = max(peak_heights);
        peak_position = peak_positions(max_idx);
        peak_height = peak_heights(max_idx);
    end
end



%% Plotter
lw = 2;
fs = 12;
alw= 1.5;
target_x = -0.33; % Target x value for peak
tolerance = 0.1; % Tolerance range around target x value
N = 10;

gr = [147, 149, 151];

figure(1)

tiledlayout(2, 1, 'TileSpacing', 'Tight', 'Padding', 'Tight');
nexttile
% Plot the bar chart
b = bar(DCF_avg, 'FaceColor', 'flat', 'FaceAlpha', 0.8); hold on 
b.CData(1,:) = [211 211 211]./255;
b.CData(2,:) = [211 211 211]./255;
b.CData(3,:) = [211 211 211]./255;
b.CData(4,:) = [211 211 211]./255;
b.CData(5,:) = [211 211 211]./255;
b.CData(6,:) = [211 211 211]./255;
b.CData(7,:) = [211 211 211]./255;
b.CData(8,:) = [211 211 211]./255;
b.CData(9,:) = [211 211 211]./255;
b.CData(10,:) = [211 211 211]./255;
errorbar(1:10, DCF_avg, DCF_std, 'k', 'linestyle', 'none', 'CapSize', 10);
hold off;

title('Intracellular [ROS]');
xlabel('Strain');
ylabel('DCF fluorescence per cell (A.U./OD_{750})');
ax = gca; 
ax.FontSize = fs;
xticks(1:10);
ax.TickLabelInterpreter = 'none';
ax.XTickLabel = titles;
ax.LineWidth = alw;
% Make edge line thicker for specific strains
highlightStrains = [1, 6];
for i = highlightStrains
    % Add a horizontal line
    x0 = i; % x value of the bar
    y0 = DCF_avg(i); % y value of the bar
    line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
end
box off


nexttile
b = bar(SW_diff, 'FaceColor', 'flat', 'FaceAlpha', 0.8);hold on
b.CData(1,:) = 'k';
b.CData(2,:) = cmapT(2,:);
b.CData(3,:) = cmapT(2,:);
b.CData(4,:) = cmapT(2,:);
b.CData(5,:) = cmapT(2,:);
b.CData(6,:) = 'k';
b.CData(7,:) = cmapT(2,:);
b.CData(8,:) = cmapT(2,:);
b.CData(9,:) = cmapT(2,:);
b.CData(10,:) = cmapT(2,:);
title('Extracellular [MV]');
xlabel('Strain');
ylabel('Fraction of total [MV] in supernatant (%)');
ax = gca; 
xticks(1:10);
ax.TickLabelInterpreter = 'none';
ax.XTickLabel = titles;
ax.FontSize = fs;
ax.LineWidth = alw;
ax.YLim = [20 80];

% Make edge line thicker for wt strains
for i = highlightStrains
    % Add a horizontal line
    x0 = i; % x value of the bar
    y0 = SW_diff(i); % y value of the bar
    line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
end
box off
