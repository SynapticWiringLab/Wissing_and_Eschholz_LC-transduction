% This script plots the efficacy of transgene expression using different model systems for LC targetting 
% as reported in Wissing and Eschholz et al, PLOS Biology, 2025, Figure 1, panel D. 
%(C)opyright Alexander Dieter, 2025


% prepare workspace
clear all; close all; clc; 

% load data (output of script "a01_quantify_efficacy_specificity_single"), mdefine experimental groups and plot colors for each group
files = {'data_DBH_50_WissingEschholz2025.mat', 'data_NET_50_WissingEschholz2025.mat', 'data_TH_50_WissingEschholz2025.mat', 'data_PRS_50_WissingEschholz2025.mat', };
grouplab = {'DBH-cre', 'NET-cre', 'TH-cre', 'PRSx8', };
groupcol = [225 031 028; 146 39 143; 244 126 032; 056 126 184] ./ 255;



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
        n_TH =          length(cell2mat(d.co_express(idx)));    % get number of TH positive neurons detected in this mouse
        n_co_THbased =  sum(cell2mat(d.co_express(idx)));       % get number of co-expressing neurons based on TH staining in this mouse
        efficacy = n_co_THbased/n_TH;                           % compute fraction of co-expressing neurons based on TH staining
        
        % Put in the table
        data(mouseIdx, groupIdx) = efficacy.*100;
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

bar(3, mean(data(:, 3)), 'FaceColor', groupcol(3, :), 'FaceAlpha', 0.5)
plot([2.8:0.4/(length(data(:, 3))-1):3.2], data(:, 3), 'ko', 'color', groupcol(3, :))
errorbar(3, mean(data(:, 3)), std(data(:, 3))./sqrt(length(data(:, 3))), 'ko-', 'color', groupcol(3, :), 'LineWidth', 2, 'MarkerFaceColor', groupcol(3, :), 'MarkerSize', 8);

bar(4, mean(data(:, 4)), 'FaceColor', groupcol(4, :), 'FaceAlpha', 0.5)
plot([3.8:0.4/(length(data(:, 4))-1):4.2], data(:, 4), 'ko', 'color', groupcol(4, :))
errorbar(4, mean(data(:, 4)), std(data(:, 4))./sqrt(length(data(:, 4))), 'ko-', 'color', groupcol(4, :), 'LineWidth', 2, 'MarkerFaceColor', groupcol(4, :), 'MarkerSize', 8);


set(gca, 'XTick', 1:4, 'XTickLabel', {'DBH-cre', 'NET-cre', 'TH-cre', 'PRSx8'})
xlim([0.5 4.5]); ylim([0 100]); ylabel('efficacy [%]')


[p, tbl, stats]  = anova1(data);
c = multcompare(stats)





