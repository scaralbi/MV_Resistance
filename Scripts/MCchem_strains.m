
%% Import Data
cd ../Data/Echem_MV

DPV = readmatrix('AS_SWV_DPV_0hr15hr.xlsx', 'Sheet', 'DPV');
SWV = readmatrix('AS_SWV_DPV_0hr15hr.xlsx', 'Sheet', 'SWV');




cd ../../Scripts


%% Params
cmap = brewermap(6,'Blues');
cmap_blanked = cmap(2:end, :);
cmapT = brewermap(3,'Set2');
N = 10;

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

titles = {'wt_Nixon','mvR01_Nixon', 'mvR02_Nixon', 'mvR03_Nixon', 'mvR06_Nixon', 'wt_Howe', 'mvR09_Howe', 'mvR10_Howe', 'mvR11_Howe', 'mvR12_Howe'};


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


SWV0_A = [SWV_strains(:,1), SWV_strains(:,3), SWV_strains(:,5), SWV_strains(:,7), SWV_strains(:,9), SWV_strains(:,11), SWV_strains(:,13), SWV_strains(:,15), SWV_strains(:,17), SWV_strains(:,19)];
SWV15_A = [SWV_strains(:,2), SWV_strains(:,4), SWV_strains(:,6), SWV_strains(:,8), SWV_strains(:,10), SWV_strains(:,12), SWV_strains(:,14), SWV_strains(:,16), SWV_strains(:,18), SWV_strains(:,20)];


SWV0 = [SWV_blanked(:,1), SWV_blanked(:,3), SWV_blanked(:,5), SWV_blanked(:,7), SWV_blanked(:,9), SWV_blanked(:,11), SWV_blanked(:,13), SWV_blanked(:,15), SWV_blanked(:,17), SWV_blanked(:,19)];
SWV15 = [SWV_blanked(:,2), SWV_blanked(:,4), SWV_blanked(:,6), SWV_blanked(:,8), SWV_blanked(:,10), SWV_blanked(:,12), SWV_blanked(:,14), SWV_blanked(:,16), SWV_blanked(:,18), SWV_blanked(:,20)];

DPV0 = [DPV_blanked(:,1), DPV_blanked(:,3), DPV_blanked(:,5), DPV_blanked(:,7), DPV_blanked(:,9), DPV_blanked(:,11), DPV_blanked(:,13), DPV_blanked(:,15), DPV_blanked(:,17), DPV_blanked(:,19)];
DPV15 = [DPV_blanked(:,2), DPV_blanked(:,4), DPV_blanked(:,6), DPV_blanked(:,8), DPV_blanked(:,10), DPV_blanked(:,12), DPV_blanked(:,14), DPV_blanked(:,16), DPV_blanked(:,18), DPV_blanked(:,20)];


DPV_area = DPV_blanked(k_dp(1):k_dp(end), :);
SWV_area = SWV_blanked(k_sw(1):k_sw(end), :);

SWV0_strains = SWV0_A(k_sw0(1):k_sw0(end), :);
SWV15_strains = SWV15_A(k_sw15(1):k_sw15(end), :);


SWV_blanked = SWV_strains - SWV_blank;

%calculate baseline from  v data
%m = (y2-y1)/(x1-x1); %for SWV
BlankLine_Y0 = [];
for s = 1: N
    m = (SWV0_strains(end,s)-SWV0_strains(1,s))/(SWV_area_V0(end,1)-SWV_area_V0(1,1));
    c = SWV0_strains(1,s) - m * SWV_area_V0(1,1);
    %c = y1 - m x1
    for i = 1:length(SWV_area_V0)
        BlankLine_Y0(i,s) = m * (SWV_area_V0(i,1)) + c;
    end
end

BlankLine_Y15 = [];

for s = 1: N
    m = (SWV15_strains(end,s)-SWV15_strains(1,s))/(SWV_area_V15(end,1)-SWV_area_V15(1,1));
    c = SWV15_strains(1,s) - m * SWV_area_V15(1,1);
    for i = 1:length(SWV_area_V15)
        BlankLine_Y15(i,s) = m * (SWV_area_V15(i,1)) + c;
    end
end


SWV0_area = SWV0(k_sw0(1):k_sw0(end), :);
SWV15_area = SWV15(k_sw15(1):k_sw15(end), :);


DPV0_area = DPV0(k_dp(1):k_dp(end), :);
DPV15_area = DPV15(k_dp(1):k_dp(end), :);

DPV_int = trapz(max(DPV_area,0));

SWV_int = trapz(max(SWV_area,0));

% SWV0_int = trapz(max(SWV0_area,0));
% SWV15_int = trapz(max(SWV15_area,0));

DPV0_int = trapz(max(DPV0_area,0));
DPV15_int = trapz(max(DPV15_area,0));


SWV0_baselined = SWV0_strains - BlankLine_Y0;
SWV15_baselined = SWV15_strains - BlankLine_Y15;

SWV0_int = trapz(max(SWV0_baselined,0));
SWV15_int = trapz(max(SWV15_baselined,0));

% 
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
fs = 14;
alw= 1.5;
target_x = -0.33; % Target x value for peak
tolerance = 0.1; % Tolerance range around target x value
% 
% figure(1)
% tiledlayout(2, 5, 'TileSpacing', 'Tight', 'Padding', 'Tight');
% 
% % Create  plots for the 2x5 layout
% peak_positions_DPV = cell(length(titles), 1);
% peak_heights_DPV = cell(length(titles), 1);
% for i = 1:N
%     nexttile
%     t0 = 2*i - 1;
%     t15 = 2*i;
%     plot(DPV_area_V, DPV_area(:,t0), 'LineWidth', lw, 'Color', cmapT(1,:)); hold on
%     [peak_position, peak_height] = find_peak_around_target(DPV_area_V, DPV_area(:,t0), target_x, tolerance);
%     peak_positions_DPV{i} = peak_position;
%     peak_heights_DPV{i} = peak_height;
%     % plot(peak_position, peak_height, 'x', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
%     %next time poiint
%     plot(DPV_area_V, DPV_area(:,t15), 'LineWidth', lw, 'Color', cmapT(3,:)); 
%     [peak_position, peak_height] = find_peak_around_target(DPV_area_V, DPV_area(:,t15), target_x, tolerance);
%     peak_positions_DPV{i} = peak_position;
%     peak_heights_DPV{i} = peak_height;
%     % plot(peak_position, peak_height, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
% 
%     ax = gca; 
%     ax.FontSize = fs;
%     ax.LineWidth = alw;
%     ax.YLim = [0.5 4.5];
%     xlabel('Voltage (V vs. SHE)');
%     ylabel('Current (\muA)');
%     % ax.XLim = [-0.8 0.5];
%     % ax.YLim = [-6 3];
%     title(titles{i}, 'Interpreter','none');
%     box off
%     if i == 10
%         x = [10, 10] ; % length
%         [legh, legObj, ~, ~] = legend('4 hours', '15 hours', 'box', 'off', 'Location', 'southeast');
%         hlegObj = findobj(legObj,'type','line');
%         lineW = 2 ; % Line Width 
%         set(hlegObj,'LineWidth',lineW);
%         legh.ItemTokenSize = x ;
%     end
% 
% end


%SWV raw
figure(1)
tiledlayout(2, 5, 'TileSpacing', 'Tight', 'Padding', 'Tight');
% Create  plots for the 2x5 layout
target_x = -0.35; % Target x value for peak
% tolerance = 0.1; % Tolerance range around target x value
peak_positions_SWV0S = cell(length(titles), 1);
peak_heights_SWV0S = cell(length(titles), 1);
peak_positions_SWV15S = cell(length(titles), 1);
peak_heights_SWV15S = cell(length(titles), 1);
for i = 1:N
    nexttile
    plot(SWV_area_V0, SWV0_strains(:,i), 'LineWidth', lw, 'Color', cmapT(1,:)); hold on
    [peak_position0S, peak_height0S] = find_peak_around_target(SWV_area_V0, SWV0_strains(:,i), target_x, tolerance);
    peak_positions_SWV0S{i} = peak_position0S;
    peak_heights_SWV0S{i} = peak_height0S;
    plot([peak_position0S peak_position0S], [0.5 peak_height0S], 'Color', cmapT(1,:), 'LineStyle', '-.', 'LineWidth', 1);
    plot([-0.7 peak_position0S], [peak_height0S peak_height0S], 'Color', cmapT(1,:),  'LineStyle', '-.', 'LineWidth', 1);


    %Threshold line
    plot(SWV_area_V0, BlankLine_Y0(:,N), '--', 'Color', [211 211 211]./255, 'LineWidth', 1); 

    % plot(peak_position0, peak_height0, 'x', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
    %next time poiint
    plot(SWV_area_V15, SWV15_strains(:,1), 'LineWidth', lw, 'Color', cmapT(3,:)); 
    [peak_position15S, peak_height15S] = find_peak_around_target(SWV_area_V15, SWV15_strains(:,i), target_x, tolerance);
    peak_positions_SWV15S{i} = peak_position15S;
    peak_heights_SWV15S{i} = peak_height15S;
    % plot(peak_position15, peak_height15, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
     %Threshold line
    plot(SWV_area_V15, BlankLine_Y15(:,N), '--', 'Color', [211 211 211]./255, 'LineWidth', 1); 
    plot([peak_position15S peak_position15S], [0.5 peak_height15S], 'Color', cmapT(3,:), 'LineStyle', '-.', 'LineWidth', 1);
    plot([-0.7 peak_position15S], [peak_height15S peak_height15S], 'Color', cmapT(3,:),  'LineStyle', '-.', 'LineWidth', 1);
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = alw;
    ax.YLim = [0.5 1.8];
    xlabel('Voltage (V vs. SHE)');
    ylabel('Current (\muA)');
    ax.XLim = [-0.7 0.1];
    % xticks(-0.7:0.1:0.1);
    % ax.XTickLabel = (-0.7:0.1:0.1);
    % ax.YLim = [-6 3];
    title(titles{i}, 'Interpreter','none');
    box off
    if i >= 1 
        x = [10, 10] ; % length
        [legh, legObj, ~, ~] = legend('4 hours', '', '', '', '15 hours', '', '', '', 'box', 'off', 'Location', 'southeast');
        hlegObj = findobj(legObj,'type','line');
        lineW = 2 ; % Line Width 
        set(hlegObj,'LineWidth',lineW);
        legh.ItemTokenSize = x ;
    end

end


% %SWV blanked
% figure(2)
% tiledlayout(2, 5, 'TileSpacing', 'Tight', 'Padding', 'Tight');
% % Create  plots for the 2x5 layout
% target_x = -0.35; % Target x value for peak
% % tolerance = 0.1; % Tolerance range around target x value
% peak_positions_SWV0 = cell(length(titles), 1);
% peak_heights_SWV0 = cell(length(titles), 1);
% peak_positions_SWV15 = cell(length(titles), 1);
% peak_heights_SWV15 = cell(length(titles), 1);
% for i = 1:N
%     nexttile
%     plot(SWV_area_V0, SWV0_baselined(:,i), 'LineWidth', lw, 'Color', cmapT(1,:)); hold on
%     [peak_position0, peak_height0] = find_peak_around_target(SWV_area_V0, SWV0_baselined(:,i), target_x, tolerance);
%     peak_positions_SWV0{i} = peak_position0;
%     peak_heights_SWV0{i} = peak_height0;
%     % plot(peak_position0, peak_height0, 'x', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
%     %next time poiint
%     plot(SWV_area_V15, SWV15_baselined(:,1), 'LineWidth', lw, 'Color', cmapT(3,:)); 
%     [peak_position15, peak_height15] = find_peak_around_target(SWV_area_V15, SWV15_baselined(:,i), target_x, tolerance);
%     peak_positions_SWV15{i} = peak_position15;
%     peak_heights_SWV15{i} = peak_height15;
%     % plot(peak_position15, peak_height15, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
% 
%     ax = gca; 
%     ax.FontSize = fs;
%     ax.LineWidth = alw;
%     ax.YLim = [0 1];
%     xlabel('Voltage (V vs. SHE)');
%     ylabel('Current (\muA)');
%     ax.XLim = [-0.7 0.1];
%     % ax.YLim = [-6 3];
%     title(titles{i}, 'Interpreter','none');
%     box off
%     if i == 10
%         x = [10, 10] ; % length
%         [legh, legObj, ~, ~] = legend('4 hours', '15 hours', 'box', 'off', 'Location', 'southeast');
%         hlegObj = findobj(legObj,'type','line');
%         lineW = 2 ; % Line Width 
%         set(hlegObj,'LineWidth',lineW);
%         legh.ItemTokenSize = x ;
%     end
% 
% end
% 

% titles
% SWV0_int
% SWV15_int


DPV_Area = [DPV0_int; DPV15_int];
SWV_Area = [SWV0_int; SWV15_int];


%SWV calibration from peak area from calibration experiments
% y(area) = 3.40x(mv conc) + 4.82
% x = (y - 4.82)/3.40 %when blank with 0
% y(area) = 3.40x(mv conc) + 4.82 %when blank with baseline

%SWV calibration from peak area
% y(area) = 3.37x(mv conc) + 6.98
% x = (y - 6.98) / 3.37



SW_MV0 = (SWV0_int-6.98)./3.37;
SW_MV15 = (SWV15_int-6.98)./3.37;
SWMV = [SW_MV0; SW_MV15];
SW_diff = (SW_MV15./SW_MV0).*100; %  fraction of the initial MV that remains in the supernatant after 15 hours expressed as percentage


DP_MV0 = (DPV0_int-2.13)./7.05;
DP_MV15 = (DPV15_int-2.13)./7.05;
DPMV = [DP_MV0; DP_MV15];
DP_diff = (DP_MV15./DP_MV0).*100; %  fraction of the initial MV that remains in the supernatant after 15 hours expressed as percentage


Diff = [DP_diff; SW_diff];

figure(3)
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Tight');
% nexttile
% bar_handle = bar(DPV_Area', 'grouped', 'FaceColor', 'flat', 'FaceAlpha', 0.6);
% bar_handle(1).CData = cmapT(1,:);
% bar_handle(2).CData = cmapT(3,:);
% 
% title('DPV - Charge');
% xlabel('Strain');
% ylabel('Peak area (\muC)')
% legend('t = 0 h','t = 15 h', 'box', 'off');
% ax = gca; 
% xticks(1:10);
% ax.TickLabelInterpreter = 'none';
% ax.XTickLabel = titles;
% ax.FontSize = fs;
% ax.LineWidth = alw;
% box off

nexttile
bar_handle = bar(SWV_Area', 'grouped', 'grouped', 'FaceColor', 'flat', 'FaceAlpha', 0.5, 'LineWidth',1, 'EdgeColor', 'k');
bar_handle(1).CData = cmapT(1,:);
bar_handle(2).CData = cmapT(3,:);

title('Total charge transferred from culture supernatant to electrodes after methyl viologen treatment');
xlabel('Strain');
ylabel('Peak area (\muC)')
legend('t = 0 h','t = 15 h', 'box', 'off');
ax = gca; 
xticks(1:10);
ax.TickLabelInterpreter = 'none';
ax.XTickLabel = titles;
ax.FontSize = fs;
ax.LineWidth = alw;
box off
% 
% nexttile
% bar_handle = bar(DPMV', 'grouped', 'grouped', 'FaceColor', 'flat',  'FaceAlpha', 0.6);
% bar_handle(1).CData = cmapT(1,:);
% bar_handle(2).CData = cmapT(3,:);
% 
% title('DPV - MV Concentration');
% xlabel('Strain');
% ylabel('[Methyl viologen] (\muM)');
% legend('t = 0 h','t = 15 h', 'box', 'off');
% ax = gca; 
% xticks(1:10);
% ax.TickLabelInterpreter = 'none';
% ax.XTickLabel = titles;
% ax.FontSize = fs;
% ax.LineWidth = alw;
% box off
% 
% nexttile
% bar_handle = bar(SWMV', 'grouped', 'grouped', 'FaceColor', 'flat', 'FaceAlpha', 0.6);
% bar_handle(1).CData = cmapT(1,:);
% bar_handle(2).CData = cmapT(3,:);
% 
% title('SWV - MV Concentration');
% xlabel('Strain');
% ylabel('[Methyl viologen] (\muM)');
% legend('t = 0 h','t = 15 h', 'box', 'off');
% ax = gca; 
% xticks(1:10);
% ax.TickLabelInterpreter = 'none';
% ax.XTickLabel = titles;
% ax.FontSize = fs;
% ax.LineWidth = alw;
% box off
% 
% nexttile
% b = bar(DP_diff, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
% b.CData(1,:) = 'k';
% b.CData(2,:) = cmapT(2,:);
% b.CData(3,:) = cmapT(2,:);
% b.CData(4,:) = cmapT(2,:);
% b.CData(5,:) = cmapT(2,:);
% b.CData(6,:) = 'k';
% b.CData(7,:) = cmapT(2,:);
% b.CData(8,:) = cmapT(2,:);
% b.CData(9,:) = cmapT(2,:);
% b.CData(10,:) = cmapT(2,:);
% title('Extracellular [MV] (DPV)');
% xlabel('Strain');
% ylabel('[MV] in supernatant (%)');
% ax = gca; 
% ax.FontSize = fs;
% xticks(1:10);
% ax.TickLabelInterpreter = 'none';
% ax.XTickLabel = titles;
% ax.LineWidth = alw;
% % Make edge line thicker for specific strains
% highlightStrains = [1, 6];
% for i = highlightStrains
%     % Add a horizontal line
%     x0 = i; % x value of the bar
%     y0 = DP_diff(i); % y value of the bar
%     line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
% end
% box off

nexttile
b = bar(SW_diff, 'FaceColor', 'flat', 'FaceAlpha', 0.4, 'LineStyle','-', 'LineWidth',1, 'EdgeColor', 'k');
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
title('Fraction of total methyl viologen detected in the supernatant after 15 hours post treatment');
xlabel('Strain');
ylabel('[MV]_{t=15}/[MV]_{t=15} (%)');
ax = gca; 
xticks(1:10);
ax.TickLabelInterpreter = 'none';
ax.XTickLabel = titles;
ax.FontSize = fs;
ax.LineWidth = alw;
% Make edge line thicker for wt strains
for i = highlightStrains
    % Add a horizontal line
    x0 = i; % x value of the bar
    y0 = SW_diff(i); % y value of the bar
    line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
end
box off
