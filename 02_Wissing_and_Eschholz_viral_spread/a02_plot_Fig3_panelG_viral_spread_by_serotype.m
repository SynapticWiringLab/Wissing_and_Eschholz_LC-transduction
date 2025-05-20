% This script plots the estimation of viral spread when using AAVs of different serotypes 
% as reported in Wissing and Eschholz et al, PLOS Biology, 2025, Figure 3, panel G. 
%(C)opyright Alexander Dieter, 2025

%% prepare workspace
clear all; close all; clc; 

AAV9_300nl = load('AAV9_300nl_WissingEschholz2025');    % load data of AAV2/9
AAV2_300nl = load('AAV2_300nl_WissingEschholz2025');    % load data of AAV2/2

%% define fluorescence thresholds (must be identical to previous script)
levels = 0.1:0.1:0.9; 


%% create figure and plot data
f0 = figure('color', [1 1 1]); hold on;


fill([levels, fliplr(levels)], [mean(AAV2_300nl.viral_spread)+std(AAV2_300nl.viral_spread), fliplr(mean(AAV2_300nl.viral_spread)-std(AAV2_300nl.viral_spread))], ...
        'k', 'EdgeColor', 'none', 'FaceColor',  [226 85 14]./255, 'FaceAlpha', 1/3); hold('on');
plot(levels, mean(AAV2_300nl.viral_spread), '-', 'Color',  [226 85 14]./255, 'LineWidth', 3); 
scatter(levels, mean(AAV2_300nl.viral_spread), repmat(50, length(levels), 1),  [226 85 14]./255, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

fill([levels, fliplr(levels)], [mean(AAV9_300nl.viral_spread)+std(AAV9_300nl.viral_spread), fliplr(mean(AAV9_300nl.viral_spread)-std(AAV9_300nl.viral_spread))], ...
        'k', 'EdgeColor', 'none', 'FaceColor', [11 0 155]./255, 'FaceAlpha', 1/3); hold('on');
plot(levels, mean(AAV9_300nl.viral_spread), '-', 'Color', [11 0 155]./255, 'LineWidth', 3); 
scatter(levels, mean(AAV9_300nl.viral_spread), repmat(50, length(levels), 1), [11 0 155]./255, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);


xlim([0.05 0.95]); ylabel('viral spread [mm²]'); xlabel('threshold [0-1]')
set(gca, 'XTick', levels, 'XTickLabel', levels); title('different serotypes')




%% statistics - perform a series of two-sample t-tests across different fluorescence thresholds

for statIdx = 1:length(levels)
    [h(statIdx), p(statIdx)] = ttest2(AAV9_300nl.viral_spread(:, statIdx), AAV2_300nl.viral_spread(:, statIdx));
end










