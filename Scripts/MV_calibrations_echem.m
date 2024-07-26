%% Import Data
cd ../Data/Echem_MV

DPV = readmatrix('DPV_SWV_MV_CalibrationCurve_JML.xlsx', 'Sheet', 'DPV');
SWV = readmatrix('DPV_SWV_MV_CalibrationCurve_JML.xlsx', 'Sheet', 'SWV');




cd ../../Scripts

%% Params
cmap = brewermap(6,'Blues');
cmap_blanked = cmap(2:end, :);

% Seaborn Tab10 palette
cmapw(1,:) = [72 118 177]./255;
cmapw(2,:) = [127 87 77]./255;

% Seaborn Tab10 palette
cmaps(1,:) = [227 126 35]./255;
cmaps(2,:) = [97 159 58]./255;
cmaps(3,:) = [184 43 44]./255;
cmaps(4,:) = [136 103 185]./255;
cmaps(5,:) = [202 118 190]./255;
cmaps(6,:) = [126 126 126]./255;
cmaps(7,:) = [187 188 56]./255;
cmaps(8,:) = [109 189 205]./255;

titles = {'wt_Nixon','mvR01_Nixon', 'mvR02_Nixon', 'mvR03_Nixon', 'mvR06_Nixon', 'wt_Howe', 'mvR09_Howe', 'mvR10_Howe', 'mvR11_Howe', 'mvR12_Howe'};

MV = [0 1 2.5 5 10];
MV_blanked = [1 2.5 5 10];

DPV_V = DPV(:,1)+0.2; %convert to SHE;
DPV_A = DPV(:,2:end);
DPV_blank = [DPV_A(:,1), DPV_A(:,1), DPV_A(:,1), DPV_A(:,1)];
DPV_MV = DPV_A(:, 2:end);

SWV_V = SWV(:,1)+0.2; %get volt data and convert to SHE

%calculate baseline from  v data
m = (1.17-0.543)/(0.194+0.6); %for SWV
BlankLine_Y = [];
c = 1.02;

for i = 1:length(SWV_V)
    BlankLine_Y(i) = m * (SWV_V(i)) + c;
end



SWV_A = SWV(:,2:end); %get y values (no voltage)
SWV_blank = [BlankLine_Y', BlankLine_Y', BlankLine_Y', BlankLine_Y'];
SWV_MV = SWV_A(:, 2:end);

%blank measurements
SWV_blanked = SWV_MV - SWV_blank;
DPV_blanked = DPV_MV - DPV_blank;

k_dp = find((-0.5<=DPV_V) & (DPV_V<=-0.2));
k_sw = find((-0.6<=SWV_V) & (SWV_V<=0.1));

DPV_area_V = DPV_V(k_dp(1):k_dp(end), 1);
SWV_area_V = SWV_V(k_sw(1):k_sw(end), 1);

DPV_area = DPV_blanked(k_dp(1):k_dp(end), :);
SWV_area = SWV_blanked(k_sw(1):k_sw(end), :);

DPV_int = trapz(max(DPV_area,0));

SWV_int = trapz(max(SWV_area,0));


 
%% Function to find peaks around a specified x value
function [peak_position, peak_height] = find_peak_around_target(voltage, current, target_x, tolerance)
    % Find peaks
    [peak_heights, peak_indices] = findpeaks(current, 'MinPeakProminence', 0.1); % Adjust prominence threshold as needed
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

%% Plot and find peaks

lw = 2;
alw= 1.5;
fs = 14;
target_x = -0.33; % Target x value for peak
tolerance = 0.1; % Tolerance range around target x value
% 
% % DPV
% figure(1)
% tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Tight');
% nexttile
% for i = 1: length(MV)
%     plot(DPV_V, DPV_A(:,i), 'Color', cmap(i,:), 'LineWidth', lw); hold on
% end
% xlabel("Voltage (V vs. SHE)");
% ylabel("Current (\muA)");
% legend('blank', 'MV = 0.1 \muM', 'MV = 1 \muM', 'MV = 2.5 \muM', 'MV = 5 \muM', 'MV = 10 \muM', "box", "off");
% ax = gca;
% ax.FontSize = fs;
% ax.XLim = [-0.8 0.2];
% ax.LineWidth = alw;
% title("Differential Pulse Voltammetry (DPV)");
% box off
% 
% % DPV blanked with peaks
% nexttile
% peak_positions_DPV = cell(length(MV_blanked), 1);
% peak_heights_DPV = cell(length(MV_blanked), 1);
% for i = 1: length(MV_blanked)
%     plot(DPV_area_V, DPV_area(:,i), 'Color', cmap_blanked(i,:), 'LineWidth', lw); hold on
%     [peak_position, peak_height] = find_peak_around_target(DPV_area_V, DPV_area(:,i), target_x, tolerance);
%     peak_positions_DPV{i} = peak_position;
%     peak_heights_DPV{i} = peak_height;
%     plot(peak_position, peak_height, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
% end
% xlabel("Voltage (V vs. SHE)");
% ylabel("Current (\muA)");
% legend('MV = 0.1 \muM', '', 'MV = 1 \muM', '', 'MV = 2.5 \muM', '', 'MV = 5 \muM', '', 'MV = 10 \muM', '', "box", "off");
% ax = gca;
% ax.LineWidth = alw;
% ax.FontSize = fs;
% ax.XLim = [-0.5 -0.2];
% title("Blank substracted and peaks detected");
% box off
% % Display the peak positions and heights
% disp('MV Peaks:');
% disp(peak_positions_DPV);
% disp(peak_heights_DPV);
% Ycal_DPV = cell2mat(peak_heights_DPV);
% 
% 
% % %% Linear Fit and Plotting for DPV Calibration
% nexttile
% scatter(MV_blanked, Ycal_DPV, 10, 'k', 'filled');
% hold on;
% 
% % Perform linear fit using fitlm
% lm_dpv = fitlm(MV_blanked, Ycal_DPV);
% 
% % Get the fitted values
% yfit_dpv = predict(lm_dpv, MV_blanked');
% 
% % Plot linear fit
% plot(MV_blanked, yfit_dpv, 'b-.', 'LineWidth', lw-1);
% 
% % Get the linear fit equation and R² value
% p_dpv = lm_dpv.Coefficients.Estimate;
% rsq_dpv = lm_dpv.Rsquared.Ordinary;
% adj_rsq_dpv = lm_dpv.Rsquared.Adjusted;
% 
% % Display linear fit equation and R² value
% eq_dpv = sprintf('y = %.2fx + %.2f\nR^2 = %.2f\nAdjusted R^2 = %.2f', p_dpv(2), p_dpv(1), rsq_dpv, adj_rsq_dpv);
% text(0.2 * max(MV_blanked), 0.9 * max(Ycal_DPV), eq_dpv, 'FontSize', 12, 'Color', 'k');
% 
% xlabel("[Methyl viologen] (\muM)");
% ylabel("Peak height (\muA)");
% ax = gca;
% ax.LineWidth = alw;
% ax.FontSize = fs;
% title("DPV calibration curve - peak height");
% box off
% hold off
% 
% 
% % Linear Fit and Plotting for DPV Calibration with Area
% nexttile
% scatter(MV_blanked, DPV_int, 10, 'k', 'filled');
% hold on;
% 
% % Perform linear fit using fitlm
% lm_dpv_int = fitlm(MV_blanked, DPV_int);
% 
% % Get the fitted values
% yfit_dpv_int = predict(lm_dpv_int, MV_blanked');
% 
% % Plot linear fit
% plot(MV_blanked, yfit_dpv_int, 'b-.', 'LineWidth', lw-1);
% 
% % Get the linear fit equation and R² value
% p_dpv_int = lm_dpv_int.Coefficients.Estimate;
% rsq_dpv_int = lm_dpv_int.Rsquared.Ordinary;
% adj_rsq_dpv_int = lm_dpv_int.Rsquared.Adjusted;
% 
% % Display linear fit equation and R² value
% eq_dpv_int = sprintf('y = %.2fx + %.2f\nR^2 = %.2f\nAdjusted R^2 = %.2f', p_dpv_int(2), p_dpv_int(1), rsq_dpv_int, adj_rsq_dpv_int);
% text(0.2 * max(MV_blanked), 0.9 * max(DPV_int), eq_dpv_int, 'FontSize', 12, 'Color', 'k');
% 
% xlabel("[Methyl viologen] (\muM)");
% ylabel("Peak area (\muC)");
% ax = gca;
% ax.LineWidth = alw;
% ax.FontSize = fs;
% title("DPV calibration curve - peak area");
% box off
% hold off;



%SWV raw
figure(2)
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Tight');

nexttile
for i = 1: length(MV)
    plot(SWV_V, SWV_A(:,i), 'Color', cmap(i,:), 'LineWidth', lw); hold on
end
%Threshold line
plot(SWV_V, BlankLine_Y, '--', 'Color', 'k', 'LineWidth', 1); 
xlabel("Voltage (V vs. SHE)");
ylabel("Current (\muA)");
legend('blank', 'MV = 1 \muM', 'MV = 2.5 \muM', 'MV = 5 \muM', 'MV = 10 \muM', '', "box", "off");
ax = gca;
ax.LineWidth = alw;
ax.FontSize = fs;
ax.XLim = [-0.6 -0.1];
title("Square Wave Voltammetry (SWV)");
box off




% SWV blanked with area and peaks
nexttile
peak_positions_SWV = cell(length(MV_blanked), 1);
peak_heights_SWV = cell(length(MV_blanked), 1);
for i = 1: length(MV_blanked)
    plot(SWV_area_V, SWV_area(:,i), 'Color', cmap_blanked(i,:), 'LineWidth', lw); hold on
    [peak_position, peak_height] = find_peak_around_target(SWV_area_V, SWV_area(:,i), target_x, tolerance);
    peak_positions_SWV{i} = peak_position;
    peak_heights_SWV{i} = peak_height;
    % plot(peak_position, peak_height, 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 2); % Optional: plot peak
    plot([peak_position peak_position], [0 peak_height], 'Color', cmap_blanked(i,:), 'LineStyle', '-.', 'LineWidth', 1);
    plot([-0.6 peak_position], [peak_height peak_height], 'Color', cmap_blanked(i,:),  'LineStyle', '-.', 'LineWidth', 1);

end
xlabel("Voltage (V vs. SHE)");
ylabel("Current (\muA)");
legend('MV = 1 \muM', '', 'MV = 2.5 \muM', '', 'MV = 5 \muM', '', 'MV = 10 \muM', '', "box", "off");
ax = gca;
ax.LineWidth = alw;
ax.FontSize = fs;
ax.XLim = [-0.6 -0.1];
ax.YLim = [0 1];
title("Blank substraction and peak detection");
box off
% 
% disp('MV peaks:');
% disp(peak_positions_SWV);
% disp(peak_heights_SWV);
% 
% Ycal_SWV = cell2mat(peak_heights_SWV);
% 
% % Linear Fit and Plotting for SWV Calibration
% nexttile
% scatter(MV_blanked, Ycal_SWV, 10, 'k', 'filled');
% hold on;
% 
% % Perform linear fit using fitlm
% lm_swv = fitlm(MV_blanked, Ycal_SWV);
% 
% % Get the fitted values
% yfit_swv = predict(lm_swv, MV_blanked');
% 
% % Plot linear fit
% plot(MV_blanked, yfit_swv, 'b-.', 'LineWidth', lw-1);
% 
% % Get the linear fit equation and R² value
% p_swv = lm_swv.Coefficients.Estimate;
% rsq_swv = lm_swv.Rsquared.Ordinary;
% adj_rsq_swv = lm_swv.Rsquared.Adjusted;
% 
% % Display linear fit equation and R² value
% eq_swv = sprintf('y = %.2fx + %.2f\nR^2 = %.2f\nAdjusted R^2 = %.2f', p_swv(2), p_swv(1), rsq_swv, adj_rsq_swv);
% text(0.2 * max(MV_blanked), 0.9 * max(Ycal_SWV), eq_swv, 'FontSize', 12, 'Color', 'k');
% xlabel("[Methyl viologen] (\muM)");
% ylabel("Peak height (\muA)");
% ax = gca;
% ax.LineWidth = alw;
% ax.FontSize = fs;
% title("SWV calibration curve - peak height");
% box off
% hold off;
% 
% 
% % Linear Fit and Plotting for SWV Calibration - Area
% nexttile
% scatter(MV_blanked, SWV_int, 20, 'k', 'filled');
% hold on;
% 
% % Perform linear fit using fitlm
% lm_swv_int = fitlm(MV_blanked, SWV_int);
% 
% % Get the fitted values
% yfit_swv_int = predict(lm_swv_int, MV_blanked');
% 
% % Plot linear fit
% plot(MV_blanked, yfit_swv_int, 'b-.', 'LineWidth', lw-1);
% 
% % Get the linear fit equation and R² value
% p_swv_int = lm_swv_int.Coefficients.Estimate;
% rsq_swv_int = lm_swv_int.Rsquared.Ordinary;
% adj_rsq_swv_int = lm_swv_int.Rsquared.Adjusted;
% 
% % Display linear fit equation and R² value
% eq_swv_int = sprintf('y = %.2fx + %.2f\nR^2 = %.2f\nAdjusted R^2 = %.2f', p_swv_int(2), p_swv_int(1), rsq_swv_int, adj_rsq_swv_int);
% text(0.2 * max(MV_blanked), 0.9 * max(SWV_int), eq_swv_int, 'FontSize', 12, 'Color', 'k');
% 
% xlabel("[Methyl viologen] (\muM)");
% ylabel("Peak area (\muC)");
% ax = gca;
% ax.LineWidth = alw;
% ax.FontSize = fs;
% title("Standard curve with peak area");
% box off
% hold off;
% 
% 
% 
