% This script plots the specificity of transgene expression using different promoters (hSyn vs CAG) for cre-dependent LC targetting in DBH-cre mice
% as reported in Wissing and Eschholz et al, PLOS Biology, 2025, Figure 3, panel D. 
%(C)opyright Alexander Dieter, 2025


% prepare workspace
clear all; close all; clc; 

% load data (output of script "a01_quantify_efficacy_specificity_single"), mdefine experimental groups and plot colors for each group
files = {'data_DBH_hSyn_50_WissingEschholz2025.mat', 'data_DBH_50_WissingEschholz2025.mat'};

grouplab = {'DBH-cre (hSyn)', 'DBH-cre (CAG)'};
groupcol = [255 197 11; 189 31 45]./ 255;



% Loop over groups
n_groups = numel(files);
for groupIdx = 1:n_groups
    
    load(files{groupIdx});    % Load preprocessed data
    
    % Loop over single mice
    mouseID = unique(d.ID);
    n_mice = numel(mouseID);
    for mouseIdx = 1:n_mice
        
        idx = (d.ID == mouseID(mouseIdx)); % Identify slices belonging to the currently considered mouse
        
        % Compute efficacy conditionally on red channel (TH+)
        n_GFP =          length(cell2mat(d.co_express_2(idx)));    % get number of TH positive neurons detected in this mouse
        n_co_GFPbased =  sum(cell2mat(d.co_express_2(idx)));       % get number of co-expressing neurons based on TH staining in this mouse
        specificity = n_co_GFPbased/n_GFP;                           % compute fraction of co-expressing neurons based on TH staining
        
        % Put in the table
        data(mouseIdx, groupIdx) = specificity.*100;
    end
end



% create plot 
figure('color', [1 1 1]); hold on;
bar(1, mean(data(:, 1)), 'FaceColor', groupcol(1, :), 'FaceAlpha', 0.5)
plot([0.8:0.4/(length(data(:, 1))-1):1.2], data(:, 1), 'ko', 'color', groupcol(1, :))
errorbar(1, mean(data(:, 1)), std(data(:, 1))./sqrt(length(data(:, 1))), 'ko-', 'color', groupcol(1, :), 'LineWidth', 2, 'MarkerFaceColor', groupcol(1, :), 'MarkerSize', 8);

bar(2, mean(data(:, 2)), 'FaceColor', groupcol(2, :), 'FaceAlpha', 0.5)
plot([1.8:0.4/(length(data(:, 2))-1):2.2], data(:, 2), 'ko', 'color', groupcol(2, :))
errorbar(2, mean(data(:, 2)), std(data(:, 2))./sqrt(length(data(:, 2))), 'ko-', 'color', groupcol(2, :), 'LineWidth', 2, 'MarkerFaceColor', groupcol(2, :), 'MarkerSize', 8);


set(gca, 'XTick', 1:2, 'XTickLabel', {'hSyn', 'CAG'})
xlim([0.5 2.5]); ylim([0 100]); ylabel('specificity [%]')


[h, p] = ttest2(data(:, 1), data(:, 2))





