
%% Import Data
cd ../Data/Echem_MV/Voltammetry_MV


SWV = readmatrix('MV_Supernatant_Samples.csv');
PH = readmatrix('pH_Supernatant_Samples.csv');
Cal = readmatrix('MV_Supernatant_Calibration.csv');
 

cd ../../../Scripts


%% Params
cmap = brewermap(3,'Blues');
cmap_blanked = cmap(2:end, :);
cmapT = brewermap(3,'Dark2');
N = 10;


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


SWV_V = SWV(4:end,1)+0.2; %convert to SHE
SWV_A = SWV(4:end,2:end);
SWV_strains = SWV_A(:, 2:end);
SWV_blank = repmat(SWV_A(:,1), 1, size(SWV_strains,2));



%blank measurements
SWV_blanked = SWV_strains - SWV_blank;

% k_sw = find((-0.5<=SWV_V) & (SWV_V<=-0.2));
k_sw = find((-0.67<=SWV_V) & (SWV_V<=-0.03));

k_sw15 = find((-0.67<=SWV_V) & (SWV_V<=-0.03));


SWV_area_V = SWV_V(k_sw(1):k_sw(end), 1);


SWV_area_V15 = SWV_V(k_sw15(1):k_sw15(end), 1);


SWV_1 = [SWV_strains(:,1), SWV_strains(:,4), SWV_strains(:,7), SWV_strains(:,10), SWV_strains(:,13), SWV_strains(:,16), SWV_strains(:,19), SWV_strains(:,22), SWV_strains(:,25), SWV_strains(:,28)];
SWV_2 = [SWV_strains(:,2), SWV_strains(:,5), SWV_strains(:,8), SWV_strains(:,11), SWV_strains(:,14), SWV_strains(:,17), SWV_strains(:,20), SWV_strains(:,23), SWV_strains(:,26), SWV_strains(:,29)];
SWV_3 = [SWV_strains(:,3), SWV_strains(:,6), SWV_strains(:,9), SWV_strains(:,12), SWV_strains(:,15), SWV_strains(:,18), SWV_strains(:,21), SWV_strains(:,24), SWV_strains(:,27), SWV_strains(:,30)];


SWV1_B = [SWV_blanked(:,1), SWV_blanked(:,4), SWV_blanked(:,7), SWV_blanked(:,10), SWV_blanked(:,13), SWV_blanked(:,16), SWV_blanked(:,19), SWV_blanked(:,22), SWV_blanked(:,25), SWV_blanked(:,28)];
SWV2_B = [SWV_blanked(:,2), SWV_blanked(:,5), SWV_blanked(:,8), SWV_blanked(:,11), SWV_blanked(:,14), SWV_blanked(:,17), SWV_blanked(:,20), SWV_blanked(:,23), SWV_blanked(:,26), SWV_blanked(:,29)];
SWV3_B = [SWV_blanked(:,3), SWV_blanked(:,6), SWV_blanked(:,9), SWV_blanked(:,12), SWV_blanked(:,15), SWV_blanked(:,18), SWV_blanked(:,21), SWV_blanked(:,24), SWV_blanked(:,27), SWV_blanked(:,30)];



SWV_area = SWV_blanked(k_sw(1):k_sw(end), :);

SWV1_strains = SWV_1(k_sw15(1):k_sw15(end), :);
SWV2_strains = SWV_2(k_sw15(1):k_sw15(end), :);
SWV3_strains = SWV_3(k_sw15(1):k_sw15(end), :);


SWV_blanked = SWV_strains - SWV_blank;

%calculate baseline from  v data
%m = (y2-y1)/(x1-x1); %for SWV
BlankLine_Y1 = [];
for s = 1: N
    m = (SWV1_strains(end,s)-SWV1_strains(1,s))/(SWV_area_V15(end,1)-SWV_area_V15(1,1));
    c = SWV1_strains(1,s) - m * SWV_area_V15(1,1);
    %c = y1 - m x1
    for i = 1:length(SWV_area_V15)
        BlankLine_Y1(i,s) = m * (SWV_area_V15(i,1)) + c;
    end
end

BlankLine_Y2 = [];
for s = 1: N
    m = (SWV2_strains(end,s)-SWV2_strains(1,s))/(SWV_area_V15(end,1)-SWV_area_V15(1,1));
    c = SWV2_strains(1,s) - m * SWV_area_V15(1,1);
    for i = 1:length(SWV_area_V15)
        BlankLine_Y2(i,s) = m * (SWV_area_V15(i,1)) + c;
    end
end

BlankLine_Y3 = [];
for s = 1: N
    m = (SWV3_strains(end,s)-SWV3_strains(1,s))/(SWV_area_V15(end,1)-SWV_area_V15(1,1));
    c = SWV3_strains(1,s) - m * SWV_area_V15(1,1);
    for i = 1:length(SWV_area_V15)
        BlankLine_Y3(i,s) = m * (SWV_area_V15(i,1)) + c;
    end
end


SWV1_area = SWV1_B(k_sw15(1):k_sw15(end), :);
SWV2_area = SWV2_B(k_sw15(1):k_sw15(end), :);
SWV3_area = SWV3_B(k_sw15(1):k_sw15(end), :);


SWV_int = trapz(max(SWV_area,0));

% SWV0_int = trapz(max(SWV0_area,0));



SWV1_baselined = SWV1_strains - BlankLine_Y1;
SWV2_baselined = SWV2_strains - BlankLine_Y2;
SWV3_baselined = SWV3_strains - BlankLine_Y3;

SWV1_int = trapz(max(SWV1_baselined,0));
SWV2_int = trapz(max(SWV2_baselined,0));
SWV3_int = trapz(max(SWV3_baselined,0));



Charges = [SWV1_int; SWV2_int; SWV3_int];
C_avg = mean(Charges,1);
C_std = std(Charges,0,1);


% pH Data
pH = PH(2:end, 2:end)';
pH_avg = mean(pH,1);
pH_std = std(pH,0,1);


%Calib data
Cal_V = Cal(4:end,1)+0.2; %convert to SHE
Cal = Cal(4:end,2:end);

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


%SWV raw
figure(1)
letters = char(65:64+(2+N));  % Where N is the number of strains

tiledlayout(3, 10, 'TileSpacing', 'Tight', 'Padding', 'compact');
nexttile([1,5])
plot(Cal_V, Cal(:,1), 'LineWidth', lw, 'Color',cmap(1,:)); hold on
plot(Cal_V, Cal(:,3), 'LineWidth', lw, 'Color',cmap(2,:)); 
plot(Cal_V, Cal(:,5), 'LineWidth', lw, 'Color',cmap(3,:)); 
legend('BG11', 'BG11 + 2 \muM MV', 'BG11 + 6 \muM MV', 'box', 'off', 'Location', 'northeast');
ax = gca; 
ax.FontSize = fs;
ax.LineWidth = alw;
% ax.YLim = [0.5 1.6];
xlabel('Voltage (V vs. SHE)');
ylabel('Current (\muA)');
addTileLabel(gca, [letters(1) ' ']);
title('Square wave voltammetry - abiotic controls')
ax.XLim = [-0.8 0.6];
% xticks(-0.7:0.1:0.1);
% ax.XTickLabel = (-0.7:0.1:0.1);
% ax.YLim = [-6 3];
box off

nexttile([1,5])
b = bar(pH_avg, 'FaceColor', 'flat', 'FaceAlpha', 0.8); hold on 
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
errorbar(1:10, pH_avg, pH_std, 'k', 'linestyle', 'none', 'CapSize', 10);
hold off;
addTileLabel(gca, [letters(2) ' ']);
title('pH of spent culture media');
xlabel('Strain');
ylabel('pH');
ax = gca; 
ax.FontSize = fs;
xticks(1:10);
ax.TickLabelInterpreter = 'none';
ax.XTickLabel = titles;
ax.LineWidth = alw;
ax.YScale = "log";
ax.YLim = [9 12];
% Make edge line thicker for specific strains
highlightStrains = [1, 6];
for i = highlightStrains
    % Add a horizontal line
    x0 = i; % x value of the bar
    y0 = pH_avg(i); % y value of the bar
    line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
end
box off



% Create  plots for the 2x5 layout
target_x = -0.35; % Target x value for peak
% tolerance = 0.1; % Tolerance range around target x value
peak_positions_SWV1S = cell(length(titles), 1);
peak_heights_SWV1S = cell(length(titles), 1);
peak_positions_SWV2S = cell(length(titles), 1);
peak_heights_SWV2S = cell(length(titles), 1);
peak_positions_SWV3S = cell(length(titles), 1);
peak_heights_SWV3S = cell(length(titles), 1);
for i = 1:N
    nexttile([1,2])
    plot(SWV_area_V15, SWV1_strains(:,i), 'LineWidth', lw, 'Color', cmapT(1,:)); hold on
    [peak_position1S, peak_height1S] = find_peak_around_target(SWV_area_V15, SWV1_strains(:,i), target_x, tolerance);
    peak_positions_SWV1S{i} = peak_position1S;
    peak_heights_SWV1S{i} = peak_height1S;
    plot([peak_position1S peak_position1S], [0.5 peak_height1S], 'Color', cmapT(1,:), 'LineStyle', '-.', 'LineWidth', 0.5);
    plot([-0.7 peak_position1S], [peak_height1S peak_height1S], 'Color', cmapT(1,:),  'LineStyle', '-.', 'LineWidth', 0.5);

    %Threshold line
    % plot(SWV_area_V15, BlankLine_Y1(:,N), '--', 'Color', [211 211 211]./255, 'LineWidth', 1);

    y1 = [SWV1_strains(1,i) SWV1_strains(end,i)];
    x1 = [SWV_area_V15(1,1) SWV_area_V15(end,1)];

    line(x1,y1, 'LineStyle', '--', 'Color', cmapT(1,:), 'LineWidth', 0.2)


    % plot(peak_position0, peak_height0, 'x', 'Color', 'b', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
    %next time poiint
    plot(SWV_area_V15, SWV2_strains(:,i), 'LineWidth', lw, 'Color', cmapT(2,:)); 
    [peak_position2S, peak_height2S] = find_peak_around_target(SWV_area_V15, SWV2_strains(:,i), target_x, tolerance);
    peak_positions_SWV2S{i} = peak_position2S;
    peak_heights_SWV2S{i} = peak_height2S;

    % plot(peak_position15, peak_height15, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
     %Threshold line
    % plot(SWV_area_V15, BlankLine_Y2(:,N), '--', 'Color', [211 211 211]./255, 'LineWidth', 1); 
    y2 = [SWV2_strains(1,i) SWV2_strains(end,i)];
    x2 = [SWV_area_V15(1,1) SWV_area_V15(end,1)];

    line(x2,y2, 'LineStyle', '--', 'Color', cmapT(2,:), 'LineWidth', 0.5)


    plot([peak_position2S peak_position2S], [0.5 peak_height2S], 'Color', cmapT(2,:), 'LineStyle', '-.', 'LineWidth', 0.5);
    plot([-0.7 peak_position2S], [peak_height2S peak_height2S], 'Color', cmapT(2,:),  'LineStyle', '-.', 'LineWidth', 0.5);


    %next time poiint
    plot(SWV_area_V15, SWV3_strains(:,i), 'LineWidth', lw, 'Color', cmapT(3,:)); 
    [peak_position3S, peak_height3S] = find_peak_around_target(SWV_area_V15, SWV3_strains(:,i), target_x, tolerance);
    peak_positions_SWV3S{i} = peak_position3S;
    peak_heights_SWV3S{i} = peak_height3S;

    % plot(peak_position15, peak_height15, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
     %Threshold line
    % plot(SWV_area_V15, BlankLine_Y3(:,N), '--', 'Color', [211 211 211]./255, 'LineWidth', 1); 
    plot([peak_position3S peak_position3S], [0.5 peak_height3S], 'Color', cmapT(3,:), 'LineStyle', '-.', 'LineWidth', 0.5);
    plot([-0.7 peak_position3S], [peak_height3S peak_height3S], 'Color', cmapT(3,:),  'LineStyle', '-.', 'LineWidth', 0.5);

    y3 = [SWV3_strains(1,i) SWV3_strains(end,i)];
    x3 = [SWV_area_V15(1,1) SWV_area_V15(end,1)];

    line(x3,y3, 'LineStyle', '--', 'Color', cmapT(3,:), 'LineWidth', 0.5)


    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = alw;
    ax.YLim = [0.5 1.6];
    xlabel('Voltage (V vs. SHE)');
    ylabel('Current (\muA)');
    ax.XLim = [-0.7 0];
    % xticks(-0.7:0.1:0.1);
    % ax.XTickLabel = (-0.7:0.1:0.1);
    % ax.YLim = [-6 3];
    addTileLabel(gca, [letters(i+2) ' ']);  % +2 because we already used A and B
    title(titles{i}, 'Interpreter','none');
    box off
    if i >= 1 
        x = [10, 10] ; % length
        [legh, legObj, ~, ~] = legend('r1', '', '', '', 'r2', '', '', '', 'r3', 'box', 'off', 'Location', 'southeast');
        hlegObj = findobj(legObj,'type','line');
        lineW = 2 ; % Line Width 
        set(hlegObj,'LineWidth',lineW);
        legh.ItemTokenSize = x ;
    end

end



% 
% 
% 
% figure(2)
% 
% tiledlayout(2, 1, 'TileSpacing', 'Tight', 'Padding', 'Tight');
% nexttile
% % Plot the bar chart
% b = bar(C_avg, 'FaceColor', 'flat', 'FaceAlpha', 0.8); hold on 
% b.CData(1,:) = [211 211 211]./255;
% b.CData(2,:) = [211 211 211]./255;
% b.CData(3,:) = [211 211 211]./255;
% b.CData(4,:) = [211 211 211]./255;
% b.CData(5,:) = [211 211 211]./255;
% b.CData(6,:) = [211 211 211]./255;
% b.CData(7,:) = [211 211 211]./255;
% b.CData(8,:) = [211 211 211]./255;
% b.CData(9,:) = [211 211 211]./255;
% b.CData(10,:) = [211 211 211]./255;
% errorbar(1:10, C_avg, C_std, 'k', 'linestyle', 'none', 'CapSize', 10);
% hold off;
% 
% title('Extracellular [MV]');
% xlabel('Strain');
% ylabel('Charge (\muC)');
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
%     y0 = C_avg(i); % y value of the bar
%     line([x0-0.3, x0 + 4.3], [y0, y0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.');
% end
% box off


function addTileLabel(tileHandle, labelText)
    % Get the position of the tile in normalized units
    pos = get(tileHandle, 'Position');
    
    % Adjust the X position to place the label outside to the left
    % Assuming 0.05 units to the left of the Y-axis
    xOffset = -0.05;  % Adjust as necessary to position outside the plot
    
    % Adjust the Y position to place the label above the X-axis units
    yOffset = pos(4) * 1.05;  % 5% above the top of the tile
    
    % Place the label
    annotation('textbox', [pos(1) + xOffset, pos(2) + yOffset, 0, 0], ...
        'String', labelText, ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 12);
end