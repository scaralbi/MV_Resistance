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