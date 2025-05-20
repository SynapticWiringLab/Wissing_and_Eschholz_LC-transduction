% This script calculates and plots the differences in normalized eGFP
% brightness between cells co-expressing TH and eGFP and cells expressing eGFP only,
% as reported in Wissing and Eschholz et al, PLOS Biology, 2025, Figure 1, panel F. 
%(C)opyright Alexander Dieter, 2025

%% clean workspace
clear all; close all; clc;

%% load data

DBHdata = load('data_DBH_50_WissingEschholz2025.mat');    % load DBH-cre data 
NETdata = load('data_NET_50_WissingEschholz2025.mat');    % load NET-cre data
THdata  = load('data_TH_50_WissingEschholz2025.mat');     % load TH-cre data
PRSdata = load('data_PRS_50_WissingEschholz2025.mat');    % load PRSx8 data


%% Sort normalized fluorescence of green channel into cells co-expressing TH and cells exclusively expressing GFP in DBH-cre mice
DBHdata.d.Green_coexpressing = [];  DBHdata.d.Green_GFPexclusive = []; 
DBHdata.d.Green_coexpressingID = [];  DBHdata.d.Green_GFPexclusiveID = []; 

for fileIdx = 1:length(DBHdata.d.ID)

    cur_green_GFPbased  =     DBHdata.d.AvBrightGreen_GFPbased{fileIdx};

    DBHdata.d.Green_coexpressing = horzcat(DBHdata.d.Green_coexpressing, cur_green_GFPbased(find(DBHdata.d.co_express_2{fileIdx} == 1)));
    DBHdata.d.Green_GFPexclusive = horzcat(DBHdata.d.Green_GFPexclusive, cur_green_GFPbased(find(DBHdata.d.co_express_2{fileIdx} == 0)));

    DBHdata.d.Green_coexpressingID = horzcat(DBHdata.d.Green_coexpressingID, repmat(DBHdata.d.ID(fileIdx), 1, sum(DBHdata.d.co_express_2{fileIdx} == 1)));
    DBHdata.d.Green_GFPexclusiveID = horzcat(DBHdata.d.Green_GFPexclusiveID, repmat(DBHdata.d.ID(fileIdx), 1, sum(DBHdata.d.co_express_2{fileIdx} == 0)));
end


%% Sort normalized fluorescence of green channel into cells co-expressing TH and cells exclusively expressing GFP in NET-cre mice
NETdata.d.Green_coexpressing = [];  NETdata.d.Green_GFPexclusive = []; 
NETdata.d.Green_coexpressingID = [];  NETdata.d.Green_GFPexclusiveID = []; 

for fileIdx = 1:length(NETdata.d.ID)

    cur_green_GFPbased  =     NETdata.d.AvBrightGreen_GFPbased{fileIdx};

    NETdata.d.Green_coexpressing = horzcat(NETdata.d.Green_coexpressing, cur_green_GFPbased(find(NETdata.d.co_express_2{fileIdx} == 1)));
    NETdata.d.Green_GFPexclusive = horzcat(NETdata.d.Green_GFPexclusive, cur_green_GFPbased(find(NETdata.d.co_express_2{fileIdx} == 0)));

    NETdata.d.Green_coexpressingID = horzcat(NETdata.d.Green_coexpressingID, repmat(NETdata.d.ID(fileIdx), 1, sum(NETdata.d.co_express_2{fileIdx} == 1)));
    NETdata.d.Green_GFPexclusiveID = horzcat(NETdata.d.Green_GFPexclusiveID, repmat(NETdata.d.ID(fileIdx), 1, sum(NETdata.d.co_express_2{fileIdx} == 0)));
end


%% Sort normalized fluorescence of green channel into cells co-expressing TH and cells exclusively expressing GFP in TH-cre mice
THdata.d.Green_coexpressing = [];  THdata.d.Green_GFPexclusive = []; 
THdata.d.Green_coexpressingID = [];  THdata.d.Green_GFPexclusiveID = []; 

for fileIdx = 1:length(THdata.d.ID)

    cur_green_GFPbased  =     THdata.d.AvBrightGreen_GFPbased{fileIdx};

    THdata.d.Green_coexpressing = horzcat(THdata.d.Green_coexpressing, cur_green_GFPbased(find(THdata.d.co_express_2{fileIdx} == 1)));
    THdata.d.Green_GFPexclusive = horzcat(THdata.d.Green_GFPexclusive, cur_green_GFPbased(find(THdata.d.co_express_2{fileIdx} == 0)));
   
    THdata.d.Green_coexpressingID = horzcat(THdata.d.Green_coexpressingID, repmat(THdata.d.ID(fileIdx), 1, sum(THdata.d.co_express_2{fileIdx} == 1)));
    THdata.d.Green_GFPexclusiveID = horzcat(THdata.d.Green_GFPexclusiveID, repmat(THdata.d.ID(fileIdx), 1, sum(THdata.d.co_express_2{fileIdx} == 0)));
end

%% Sort normalized fluorescence of green channel into cells co-expressing TH and cells exclusively expressing GFP in PRSx8 mice
PRSdata.d.Green_coexpressing = [];  PRSdata.d.Green_GFPexclusive = []; 
PRSdata.d.Green_coexpressingID = [];  PRSdata.d.Green_GFPexclusiveID = []; 

for fileIdx = 1:length(PRSdata.d.ID)

    cur_green_GFPbased  =     PRSdata.d.AvBrightGreen_GFPbased{fileIdx};

    PRSdata.d.Green_coexpressing = horzcat(PRSdata.d.Green_coexpressing, cur_green_GFPbased(find(PRSdata.d.co_express_2{fileIdx} == 1)));
    PRSdata.d.Green_GFPexclusive = horzcat(PRSdata.d.Green_GFPexclusive, cur_green_GFPbased(find(PRSdata.d.co_express_2{fileIdx} == 0)));
    
    PRSdata.d.Green_coexpressingID = horzcat(PRSdata.d.Green_coexpressingID, repmat(PRSdata.d.ID(fileIdx), 1, sum(PRSdata.d.co_express_2{fileIdx} == 1)));
    PRSdata.d.Green_GFPexclusiveID = horzcat(PRSdata.d.Green_GFPexclusiveID, repmat(PRSdata.d.ID(fileIdx), 1, sum(PRSdata.d.co_express_2{fileIdx} == 0)));

end


%% create vectors of normalized brightness values and corresponding labels
norm_brightness = [ DBHdata.d.Green_coexpressing DBHdata.d.Green_GFPexclusive, ...
                    NETdata.d.Green_coexpressing NETdata.d.Green_GFPexclusive, ...
                    THdata.d.Green_coexpressing THdata.d.Green_GFPexclusive, ...
                    PRSdata.d.Green_coexpressing PRSdata.d.Green_GFPexclusive]';

brightness_label = [repmat({'DBH co'}, length(DBHdata.d.Green_coexpressing), 1); repmat({'DBH ex'}, length(DBHdata.d.Green_GFPexclusive), 1); ...
                    repmat({'NET co'}, length(NETdata.d.Green_coexpressing), 1); repmat({'NET ex'}, length(NETdata.d.Green_GFPexclusive), 1); ...
                    repmat({'TH co'}, length(THdata.d.Green_coexpressing), 1); repmat({'TH ex'}, length(THdata.d.Green_GFPexclusive), 1); ...
                    repmat({'PRS co'}, length(PRSdata.d.Green_coexpressing), 1); repmat({'PRS ex'}, length(PRSdata.d.Green_GFPexclusive), 1)];


%% create average of relative eGFP expression per animal
DBH_IDs = unique(DBHdata.d.ID);
for DBHidx = 1:length(DBH_IDs)
   DBH_coexpress_means(DBHidx) = mean(DBHdata.d.Green_coexpressing(DBHdata.d.Green_coexpressingID == DBH_IDs(DBHidx)));
   DBH_GFPexclusive_means(DBHidx) = mean(DBHdata.d.Green_GFPexclusive(DBHdata.d.Green_GFPexclusiveID == DBH_IDs(DBHidx)));
end

NET_IDs = unique(NETdata.d.ID);
for NETidx = 1:length(NET_IDs)
   NET_coexpress_means(NETidx) = mean(NETdata.d.Green_coexpressing(NETdata.d.Green_coexpressingID == NET_IDs(NETidx)));
   NET_GFPexclusive_means(NETidx) = mean(NETdata.d.Green_GFPexclusive(NETdata.d.Green_GFPexclusiveID == NET_IDs(NETidx)));
end

TH_IDs = unique(THdata.d.ID);
for THidx = 1:length(TH_IDs)
   TH_coexpress_means(THidx) = mean(THdata.d.Green_coexpressing(THdata.d.Green_coexpressingID == TH_IDs(THidx)));
   TH_GFPexclusive_means(THidx) = mean(THdata.d.Green_GFPexclusive(THdata.d.Green_GFPexclusiveID == TH_IDs(THidx)));
end

PRS_IDs = unique(PRSdata.d.ID);
for PRSidx = 1:length(PRS_IDs)
   PRS_coexpress_means(PRSidx) = mean(PRSdata.d.Green_coexpressing(PRSdata.d.Green_coexpressingID == PRS_IDs(PRSidx)));
   PRS_GFPexclusive_means(PRSidx) = mean(PRSdata.d.Green_GFPexclusive(PRSdata.d.Green_GFPexclusiveID == PRS_IDs(PRSidx)));
end

%% calculate the relative difference of GFP fluorescence in TH-co-expressing cells vs. GFP-only cells
DBH_relFluoDiff = (DBH_coexpress_means-DBH_GFPexclusive_means);
NET_relFluoDiff = (NET_coexpress_means-NET_GFPexclusive_means);
TH_relFluoDiff  = (TH_coexpress_means-TH_GFPexclusive_means);
PRS_relFluoDiff = (PRS_coexpress_means-PRS_GFPexclusive_means);



%% create figure of relative eGFP expression based on TH co-labelling per animal
f2 = figure('color', [1 1 1]); hold on;

bar(1, mean(DBH_relFluoDiff), 'FaceColor', [189 031 045]./255)
bar(2, mean(NET_relFluoDiff), 'FaceColor', [146 039 143]./255)
bar(3, mean(TH_relFluoDiff), 'FaceColor', [249 146 022]./255)
bar(4, mean(PRS_relFluoDiff), 'FaceColor', [000 174 240]./255)

plot([0.8:0.4./(length(DBH_relFluoDiff)-1):1.2], DBH_relFluoDiff, 'ko', 'color', [189 031 045]./255)
plot([1.8:0.4./(length(NET_relFluoDiff)-1):2.2], NET_relFluoDiff, 'ko', 'color', [146 039 143]./255)
plot([2.8:0.4./(length(TH_relFluoDiff)-1):3.2], TH_relFluoDiff, 'ko', 'color', [249 146 022]./255)
plot([3.8:0.4./(length(PRS_relFluoDiff)-1):4.2], PRS_relFluoDiff, 'ko', 'color', [000 174 240]./255)
errorbar(1, mean(DBH_relFluoDiff), std(DBH_relFluoDiff), 'ko-', 'color', [189 031 045]./255, 'LineWidth', 2, 'MarkerSize', 8)
errorbar(2, mean(NET_relFluoDiff), std(NET_relFluoDiff), 'ko-', 'color', [146 039 143]./255, 'LineWidth', 2, 'MarkerSize', 8)
errorbar(3, mean(TH_relFluoDiff), std(TH_relFluoDiff), 'ko-', 'color', [249 146 022]./255, 'LineWidth', 2, 'MarkerSize', 8)
errorbar(4, mean(PRS_relFluoDiff), std(PRS_relFluoDiff), 'ko-', 'color', [000 174 240]./255, 'LineWidth', 2, 'MarkerSize', 8)

plot([0.5 4.5], [0 0], 'k--')

xlim([0.5 4.5]); ylim([-0.4 0.2])


%% do statistics on relative GFP expression based on TH co-labelling per animal
[p, tbl,stats]  = anova1([DBH_relFluoDiff' NET_relFluoDiff' TH_relFluoDiff' PRS_relFluoDiff']);
c = multcompare(stats)




