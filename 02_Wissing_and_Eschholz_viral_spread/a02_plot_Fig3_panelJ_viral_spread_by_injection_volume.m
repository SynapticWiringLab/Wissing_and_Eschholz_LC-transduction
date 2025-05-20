% This script plots the estimation of viral spread when using AAVs of different serotypes 
% as reported in Wissing and Eschholz et al, PLOS Biology, 2025, Figure 3, panel J. 
%(C)opyright Alexander Dieter, 2025

%% prepare workspace
clear all; close all; clc; 

AAV9_50nl = load('AAV9_50nl_WissingEschholz2025');      % load data of AAV2/9, 50 nl
AAV9_100nl = load('AAV9_100nl_WissingEschholz2025');    % load data of AAV2/9, 100 nl
AAV9_300nl = load('AAV9_300nl_WissingEschholz2025');    % load data of AAV2/9, 300 nl



%% define fluorescence thresholds (must be identical to previous script)
levels = 0.1:0.1:0.9; 


%% create figure and plot data
f0 = figure('color', [1 1 1]); hold on;
fill([levels, fliplr(levels)], [mean(AAV9_50nl.viral_spread)+std(AAV9_50nl.viral_spread), fliplr(mean(AAV9_50nl.viral_spread)-std(AAV9_50nl.viral_spread))], ...
        'k', 'EdgeColor', 'none', 'FaceColor', [0 191 191]./255, 'FaceAlpha', 1/5); hold('on');
plot(levels, mean(AAV9_50nl.viral_spread), '-', 'Color', [0 191 191]./255, 'LineWidth', 3); 
scatter(levels, mean(AAV9_50nl.viral_spread), repmat(50, length(levels), 1), [0 191 191]./255, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

fill([levels, fliplr(levels)], [mean(AAV9_100nl.viral_spread)+std(AAV9_100nl.viral_spread), fliplr(mean(AAV9_100nl.viral_spread)-std(AAV9_100nl.viral_spread))], ...
        'k', 'EdgeColor', 'none', 'FaceColor', [0 100 191]./255, 'FaceAlpha', 1/5); hold('on');
plot(levels, mean(AAV9_100nl.viral_spread), '-', 'Color', [0 100 191]./255, 'LineWidth', 3); 
scatter(levels, mean(AAV9_100nl.viral_spread), repmat(50, length(levels), 1), [0 100 191]./255, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

fill([levels, fliplr(levels)], [mean(AAV9_300nl.viral_spread)+std(AAV9_300nl.viral_spread), fliplr(mean(AAV9_300nl.viral_spread)-std(AAV9_300nl.viral_spread))], ...
        'k', 'EdgeColor', 'none', 'FaceColor', [11 0 155]./255, 'FaceAlpha', 1/5); hold('on');
plot(levels, mean(AAV9_300nl.viral_spread), '-', 'Color', [11 0 155]./255, 'LineWidth', 3); 
scatter(levels, mean(AAV9_300nl.viral_spread), repmat(50, length(levels), 1), [11 0 155]./255, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

xlim([0.05 0.95]); ylabel('viral spread [mm²]'); xlabel('threshold [0-1]')
set(gca, 'XTick', levels, 'XTickLabel', levels); title('different injection volumes')





%% statistics - perform a series of anovas/tukey tests across different fluorescence thresholds

for statIdx = 1:length(levels)

    anova_vals = [AAV9_300nl.viral_spread(:, statIdx ); AAV9_100nl.viral_spread(:, statIdx ); AAV9_50nl.viral_spread(:, statIdx )];
    anova_labels = [ones(size(AAV9_300nl.viral_spread, 1), 1); ones(size(AAV9_100nl.viral_spread, 1), 1).*2; ones(size(AAV9_50nl.viral_spread, 1), 1).*3];

    [p, tbl, stats] = anova1(anova_vals, anova_labels);
    anova_stats = multcompare(stats);
    all_pvals(1, statIdx) = anova_stats(1, end);
    all_pvals(2, statIdx) = anova_stats(2, end);
    all_pvals(3, statIdx) = anova_stats(3, end);

end









