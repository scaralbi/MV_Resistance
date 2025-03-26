%% Import Data
cd ../Data/Echem_MV

CA_BG = readmatrix('chronoamp_MVadapt_03V_multiuemstat_blockBPV_150uE_2hdl.xlsx', 'Sheet', 'BG11');
CA_MV = readmatrix('chronoamp_MVadapt_03V_multiuemstat_blockBPV_150uE_2hdl.xlsx', 'Sheet', '6uM_MV');

cd ../../Scripts

%% Params

cmap = brewermap(4,'RdYlBu');
cmaps = brewermap(8,'RdYlBu');
cmapw = brewermap(3,'RdYlBu');

%Seaborn Tab10 palette
cmapw(1,:) = [72 118 177]./255;
cmapw(2,:) = [127 87 77]./255;

%Seaborn Tab10 palette
cmaps(1,:) = [227 126 35]./255;
cmaps(2,:) = [97 159 58]./255;
cmaps(3,:) = [184 43 44]./255;
cmaps(4,:) = [136 103 185]./255;
cmaps(5,:) = [202 118 190]./255;
cmaps(6,:) = [126 126 126]./255;
cmaps(7,:) = [187 188 56]./255;
cmaps(8,:) = [109 189 205]./255;

cmap(1,:) = [227 126 35]./255;
cmap(2,:) = [227 126 35]./255;
cmap(3,:) = [136 103 185]./255;
cmap(4,:) = [136 103 185]./255;

titles = {'WT "Nixon"','mvR1', 'mvR2', 'mvR3', 'mvR6', 'WT "Howe"', 'mvR9', 'mvR10', 'mvR11', 'mvR12'};


%% Dark light panels with different times

light = 2;%duration of light period in hours 
dark = 2;%duration of dark period in hours 

max_time = 50;
number_of_patches = ceil(max_time / (dark+light) *2);


X = [];
for i = 1:number_of_patches
    if i ==1 
        X(:,1) = [0; light; light; 0];
    elseif mod(i,2)==1 % isodd
        X(:,i) = X(:,i-1) + [dark; light; light; dark];%my cycle started always with dark adaptation at start
    elseif mod(i,2)==0 % iseven
        X(:,i) = X(:,i-1) + [light; dark; dark; light];
    end
end

Y = [];
for i = 1:number_of_patches
    Y(:,i) = [0.01; 0.01; 20; 20];
end


C = zeros(1,number_of_patches,3);
for i = 1:number_of_patches
    if mod(i,2)==1 % isodd 
        C(1,i,:) = 0; %white for light
    elseif mod(i,2)==0 % iseven
        C(1,i,:) = 1; %black for dark
    end
end

%% Plotter 


lw = 2;
alw = 2;
fs = 12;
letters = char(65:64+(N));  % Where N is the number of strains

% Create a new figure
hFig = figure;
figure(1)
% Create a 2x1 tiled layout to separate the top 2x2 and bottom 2x5 grids
t1 = tiledlayout(5, 2, 'TileSpacing', 'Tight', 'Padding', 'Tight');
% Create  plots for the 2x2 layout
%Plot 1,1
nexttile
p1 = patch(X,Y,C,'FaceAlpha', 0.17, 'EdgeAlpha', 0);hold on
plot(CA_BG(:,1),CA_BG(:,2:4), 'LineWidth', lw, 'Color', cmapw(1,:));
plot(CA_MV(:,1),CA_MV(:,2:4), 'LineWidth', lw, 'Color', cmaps(2,:));
xlabel('time (hours)');
ylabel('current (\muA)');
title('wt_Nixon', 'Interpreter','none');
ax = gca;
ax.XLim = [0 21];
legend('','-MV', '', '', '+MV', '', '', 'Box', 'off');
ax.YLim = [0.05 12];
ax.XTick = 0:2:22;
ax.YTick = [0.1 1 10];
ax.YTickLabel = [0.1 1 10];
ax.XTickLabel = {0, '', 4, '', 8, '', 12, '', 16,  '', 20};
ax.YScale = "log";
addTileLabel(gca, [letters(1) ' ']);
box off
ax.LineWidth = alw;
ax.FontSize = fs;

nexttile
p1 = patch(X,Y,C,'FaceAlpha', 0.17, 'EdgeAlpha', 0);hold on
plot(CA_BG(:,13),CA_BG(:,14:16), 'LineWidth', lw, 'Color', cmapw(1,:));
plot(CA_MV(:,13),CA_MV(:,14:16), 'LineWidth', lw, 'Color', cmaps(2,:));
xlabel('time (hours)');
ylabel('current (\muA)');
title('wt_Howe','Interpreter','none');
ax = gca;
ax.XLim = [0 21];
legend('', '-MV', '', '', '+MV', '', '', 'Box', 'off');
ax.YLim = [0.05 12];
ax.XTick = 0:2:22;
ax.YTick = [0.1 1 10];
ax.YTickLabel = [0.1 1 10];
ax.XTickLabel = {0, '', 4, '', 8, '', 12, '', 16,  '', 20};
ax.YScale = "log";
addTileLabel(gca, [letters(2) ' ']);

box off
ax.LineWidth = alw;
ax.FontSize = fs;

%MVR Strains
Yindexes = [6 18 8 20 10 22 12 24];
Xindexes = [5 17 7 19 9 21 11 23];


% Yindexes = [6 8 10 12 18 20 22 24];
% Xindexes = [5 7 9 11 17 19 21 23];
MV_Strains = {'mvR01_Nixon', 'mvR02_Nixon', 'mvR03_Nixon', 'mvR06_Nixon', 'mvR09_Howe', 'mvR10_Howe', 'mvR11_Howe', 'mvR12_Howe'};
for i = 1:N-2
    nexttile
    p1 = patch(X,Y,C,'FaceAlpha', 0.17, 'EdgeAlpha', 0); hold on
    plot(CA_BG(:,Xindexes(i)),CA_BG(:,Yindexes(i)), 'LineWidth', lw, 'Color', cmapw(1,:));
    plot(CA_MV(:,Xindexes(i)),CA_MV(:,Yindexes(i)), 'LineWidth', lw, 'Color', cmaps(2,:)); 
    xlabel('time (hours)');
    ylabel('current (\muA)');
    title(MV_Strains{i}, 'Interpreter','none');
    ax = gca;
    ax.XLim = [0 21];
    legend('','-MV', '+MV', 'Box', 'off');
    ax.YLim = [0.05 12];
    ax.XTick = 0:2:22;
    ax.YTick = [0.1 1 10];
    ax.YTickLabel = [0.1 1 10];
    ax.XTickLabel = {0, '', 4, '', 8, '', 12, '', 16,  '', 20};
    ax.YScale = "log";
    addTileLabel(gca, [letters(2+i) ' ']);
    box off
    ax.LineWidth = alw;
    ax.FontSize = fs;
end
