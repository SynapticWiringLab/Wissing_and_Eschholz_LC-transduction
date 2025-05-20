% This script investigates the effect of staining on the number of
% detected cells and their fluorescence.
% 
% Copyright (c) Maxime Maheu, 2023

%% LOAD AND ANALYSE DATA
%  =====================

% Move to the home directory
workingpath = 'G:\--- common\manuscripts\Wissing_et_al_LCTransduction\PLOS_Biology\03_revision_02\code\Figure_S4'; % change to path in which data (i.e. the output of script "a01_quantify_efficacy_specificity_single_for_staining_effects" is stored)
cd(workingpath);

% Specify groups
files = dir(fullfile(workingpath, '*data*.mat'));
files = {files.name};
grouplab = {'DBH-cre', 'NET-cre', 'PRSx8', 'TH-cre'};
groupcol = [225 031 028; 146 39 143; 056 126 184; 244 126 032] ./ 255;

% Define function to fit to fluorescence distribution
fun = @(b,g,x) normcdf(b + g .* norminv(x));
Niter = 10;

% Pre-allocate output variables
Ng = numel(files);
nGFPst   = cell(1,Ng);
nGFPnt   = cell(1,Ng);
colocgnt = cell(1,Ng);
colocgst = cell(1,Ng);
fluost   = cell(1,Ng);
fluont   = cell(1,Ng);
params   = cell(1,Ng);

% Loop over groups
for ig = 1:Ng
    
    % Load preprocessed data
    load(fullfile(workingpath, files{ig}));
    
    % Normalise cell fluorescence within each slice to account for
    % difference in acquisition parameters
    norm = @(x) (x - min(x)) ./ (max(x) - min(x));
    normnt = cellfun(norm, d.AvBrightGreen_GFPbased, 'uni', 0);
    normst = cellfun(norm, d.AvBrightRed_GFPbased, 'uni', 0);
    
    % Loop over single mice
    mouseID = unique(d.ID);
    Nmice = numel(mouseID);
    for im = 1:Nmice
        
        % Identify LCs belonging to the currently considered mouse
        idx = (d.ID == mouseID(im));
        
        % Get number of GFP native cells vs. stained cells
        nGFPnt{ig}(im) = mean(cellfun(@numel, d.co_express_2(idx))); % green
        nGFPst{ig}(im) = mean(cellfun(@numel, d.co_express  (idx))); % red
        
        % Compute cell co-localization
        % - based on native eGFP segmentation
        % - based on stained GFP segmentation
        n = cell2mat(d.co_express_2(idx));
        colocgnt{ig}(im) = sum(n) ./ numel(n) .* 100;
        n = cell2mat(d.co_express(idx));
        colocgst{ig}(im) = sum(n) ./ numel(n) .* 100;
        
        % Get normalised fluorescence
        % N.B. #1 we take conditionnally on the red channel which should
        % contains more cells than the green channel
        % N.B. #2 we take only cells that are detected in both channels to
        % avoid that differences in fluorescence are caused by differences
        % in cell counts
        in = {d.co_express_2(idx), 'uni', 0};
        fluont{ig}{im} = cell2mat(cellfun(@(x,i) x(logical(i)), normnt(idx), in{:}));
        fluost{ig}{im} = cell2mat(cellfun(@(x,i) x(logical(i)), normst(idx), in{:}));
        
        % Fit function to fluorescence levels
        params{ig}(im,:) = fit_function(fun, fluont{ig}{im}, fluost{ig}{im}, Niter);
    end
end

% Get total number of mice across all groups
Ntotalmice = sum(cellfun(@numel, nGFPst));

% Average cell counts over mice
avgnGFPst = mean(cell2mat(nGFPst));
semnGFPst = std(cell2mat(nGFPst)) ./ sqrt(Ntotalmice);
avgnGFPnt = mean(cell2mat(nGFPnt));
semnGFPnt = std(cell2mat(nGFPnt)) ./ sqrt(Ntotalmice);

% Compute difference in cell counts
cellcountdif = cell2mat(nGFPst) - cell2mat(nGFPnt);

% Average co-localisation
avgcolocnt = mean(cell2mat(colocgnt));
semcolocnt = std(cell2mat(colocgnt)) ./ sqrt(Ntotalmice);
avgcolocst = mean(cell2mat(colocgst));
semcolocst = std(cell2mat(colocgst)) ./ sqrt(Ntotalmice);

% Compute difference in co-localisation
colocdif = cell2mat(colocgst) - cell2mat(colocgnt);

% Average parameter value over mice
allparams = cell2mat(params(:));
avgparams = mean(allparams, 1);
semparams = std(allparams, [], 1) ./ sqrt(Ntotalmice);

% Compute 2D distribution of fluorescence in native/stained channels
nbins = 50;
edges = linspace(0, 1, nbins);
allfluont = cell2mat(cellfun(@cell2mat, fluont, 'uni', 0));
allfluost = cell2mat(cellfun(@cell2mat, fluost, 'uni', 0));
[map, ctrs] = hist3([allfluont(:) allfluost(:)], 'Edges', {edges, edges});
k = 10;
K = (1/(k^2)) * ones(k);
map = conv2(map, K, 'same')';
map = map ./ sum(map(:));

% Compute distribution of difference between fluorescence
fluodif = allfluost - allfluont;


% % Compute distribution of difference between fluorescence per animal - added by AD, 31.10.2023
% fluodif = allfluost - allfluont;
% 
% % Loop over groups
% for ig = 1:Ng
%     
%     % Loop over single mice
%     Nmice = numel(fluost{ig});
% 
%     for im = 1:Nmice
% 
%         fluodiff{ig}{im} = fluost{ig}{im} - fluont{ig}{im};
% 
%         meanfluost{ig}(im) = mean(fluost{ig}{im});
%         meanfluont{ig}(im) = mean(fluont{ig}{im});
% 
%     end
% 
% end
% 


%% FIGURE A: CELL COUNTS IN GFP NATIVE VS. STAINED CELLS
%  =====================================================

% Difference in number of native-eGFP vs. stained eGFP cells
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare new figure
figure('Units', 'Pixels', 'Position', [77 486 701 431]);
tiledlayout(3,2);
nexttile([1,1]); hold('on');

% Plot distribution of difference between cell counts in the two channels
hx = -50:10:50;
hy = hist(cellcountdif, hx);
bar(hx, hy, 1, 'FaceColor', ones(1,3)./2, 'EdgeColor', 'None');

% Plot kernel density
ksx = -50:0.01:50;
ksy = ksdensity(cellcountdif, ksx);
plot(ksx, ksy ./ max(ksy) .* max(hy), 'k-', 'LineWidth', 1.5);

% Plot null effect
plot(zeros(1,2), [0,max(hy)], 'k:');

% Customize plot
set(gca, 'XDir', 'Reverse', 'XLim', ksx([1,end]), ...
    'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
xlabel('difference');
ylabel('counts');

% Number of native-eGFP vs. stained eGFP cells
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plot diagonal
nexttile(3,[2,1]); hold('on');
plot([0,200], [0,200], 'k--');

% Plot single mice
lgd = NaN(1,Ng);
for ig = 1:Ng
    lgd(ig) = scatter(nGFPnt{ig}, nGFPst{ig}, 50, 'filled', ...
        'MarkerFaceColor', groupcol(ig,:), 'MarkerFaceAlpha', 0.4);
end

% Plot group average
plot(avgnGFPnt + semnGFPnt .* [-1,1], repmat(avgnGFPst, 1, 2), 'k-');
plot(repmat(avgnGFPnt, 1, 2), avgnGFPst + semnGFPst .* [-1,1], 'k-');
plot(avgnGFPnt, avgnGFPst, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', ones(1,3)./2);

% Customize the plot
axis('square', 'equal');
axis(repmat([0,200], 1, 2));
set(gca, 'XTick', 0:50:200, 'YTick', 0:50:200, ...
    'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
legend(lgd, grouplab, 'Location', 'SouthEast');
xlabel('number of native eGFP^{+} cells');
ylabel('number of GFP-stained^{+} cells');

% Difference in co-localisation conditional on native-eGFP vs. stained eGFP channels
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plot distribution of difference between cell counts in the two channels
nexttile([1,1]); hold('on');
hx = -50:5:50;
hy = hist(colocdif, hx);
bar(hx, hy, 1, 'FaceColor', ones(1,3)./2, 'EdgeColor', 'None');

% Plot kernel density
ksx = -50:0.01:50;
ksy = ksdensity(colocdif, ksx);
plot(ksx, ksy ./ max(ksy) .* max(hy), 'k-', 'LineWidth', 1.5);

% Plot null effect
plot(zeros(1,2), [0,max(hy)], 'k:');

% Customize plot
set(gca, 'XDir', 'Reverse', 'XLim', ksx([1,end]), ...
    'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
xlabel('difference');
ylabel('counts');

% Co-localisation of native-eGFP vs. stained eGFP cells
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plot diagonal
nexttile([2,1]); hold('on');
plot([0,450], [0,450], 'k--');

% Plot individual mice
for ig = 1:Ng
    scatter(colocgnt{ig}, colocgst{ig}, 50, 'filled', ...
        'MarkerFaceColor', groupcol(ig,:), 'MarkerFaceAlpha', 0.4);
end

% Plot group average
plot(avgcolocnt + semcolocnt .* [-1,1], repmat(avgcolocst, 1, 2), 'k-');
plot(repmat(avgcolocnt, 1, 2), avgcolocst + semcolocst .* [-1,1], 'k-');
plot(avgcolocnt, avgcolocst, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', ones(1,3)./2);

% Customize the plot
axis('square');
set(gca, 'XTick', 0:20:100, 'YTick', 0:20:100, 'XLim', [0,100], ...
    'YLim', [0,100], 'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
xlabel('co-localization [% of native eGFP^{+}]');
ylabel('co-localization [% of GFP-stained^{+}]');

%% FIGURE B: FLUORESCENCE IN NATIVE VS. STAINED CELLS
%  ==================================================

% Difference in fluorescence in native-eGFP vs. stained eGFP channels
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare new figure
figure('Units', 'Pixels', 'Position', [779 486 528 431]);
tiledlayout(3,3);
nexttile([1,2]); hold('on');

% Plot distribution of difference between cell counts in the two channels
hx = -0.5:0.05:0.5;
hy = hist(fluodif, hx);
bar(hx, hy, 1, 'FaceColor', ones(1,3)./2, 'EdgeColor', 'None');

% Plot kernel density
ksx = -0.6:0.01:0.6;
ksy = ksdensity(fluodif, ksx);
plot(ksx, ksy ./ max(ksy) .* max(hy), 'k-', 'LineWidth', 1.5);

% Plot null effect
plot(zeros(1,2), [0,max(hy)], 'k:');

% Customize plot
set(gca, 'XDir', 'Reverse', 'XLim', ksx([1,end]), ...
    'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
xlabel('difference');
ylabel('counts');

% Distribution of fluorescence in native-eGFP vs. stained eGFP channels
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plot distribution
nexttile([2,2]); hold('on');
imagesc(edges, edges, map);
contour(edges, edges, map, 'k-');

% Plot diagonal
plot([0,1], [0,1], 'k:');

% Customize the plot
cbr = colorbar('Location', 'North');
cbr.Label.String = 'frequency [%]';
set(gca, 'XTick', 0:0.25:1, 'YTick', 0:0.25:1, ...
    'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
axis([0,1,0,1]); axis('square', 'xy');
xlabel('native eGFP fluorescence [%]');
ylabel('GFP-stained fluorescence [%]');

% Parameters describing the relation between fluorescence in both channels
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For each parameter
xjit = linspace(-0.25, 0.25, Ng);
nexttile(6,[2,1]); hold('on');
for ip = 1:2
    
    % Plot null value
    plot(0.5+[-1,0]+ip, repmat(ip-1,1,2), 'k-');
    
    % Plot individual mice
    for ig = 1:Ng
        scatter(ip+xjit(ig), params{ig}(:,ip), 50, 'filled', ...
            'MarkerFaceColor', groupcol(ig,:), 'MarkerFaceAlpha', 0.4);
    end
    
    % Plot group average
    plot(repmat(ip,1,2), avgparams(ip)+semparams(ip).*[-1,1], 'k-');
    plot(ip, avgparams(ip), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', ones(1,3)./2);
end

% Customize the plot
set(gca, 'XLim', [0,3], 'YLim', [-0.5,1.5], 'XTick', 1:2, 'XTickLabel', ...
    {'bias', 'gain'}, 'FontSize', 15, 'TickDir', 'Out', 'Box', 'Off');
ylabel('parameter value');

%% FIGURE C: EXAMPLE PARAMETER VARIATIONS
%  ======================================

% Define parameter variations
hy = {'bias variation', 'gain variation'};
gain = {[1 1 1 1 1],             [0.6 0.8 1.0 1.25 1.6]};
bias = {[-0.5 -0.25 0 0.25 0.5], [1 1 1 1 1]};
Nparam = numel(gain{1});

% Prepare figure
cmap = linspace(0, 0.75, ceil(Nparam/2))';
cmap = repmat(cat(1, flipud(cmap), cmap(2:end)), 1, 3);
xgrid = 0:0.01:1;
figure('Units', 'Pixels', 'Position', [1308 486 272 431]);

% Loop over variations
for isp = 1:2
    subplot(2,1,isp);
    
    % Display diagonals
    plot([0,1], [0,1], 'k--'); hold('on');
    plot([0,1], [1,0], 'k--');
    
    % Loop over parameter grid
    for ip = 1:Nparam
        
        % Display model predictions
        y = fun(bias{isp}(ip), gain{isp}(ip), xgrid);
        plot(xgrid, y, '-', 'Color', cmap(ip,:), 'LineWidth', 3);
    end
    
    % Customize plot
    axis('square');
    set(gca, 'XTick', [0,1], 'YTick', [0,1], 'FontSize', 15, ...
        'TickDir', 'Out', 'Box', 'Off');
    xlabel('native eGFP [%]');
    ylabel('stained GFP [%]');
    title(hy{isp}, 'FontWeight', 'Normal');
end

%% SUBFUNCTION
%  ===========

function [params, fitout, gof] = fit_function(fun, x, y, Niter)

% Create fit object from provided objective function
fitfun = fittype(fun);

% Run several iterations of the fitting procedure, varying starting
% points
fitout = cell(1,Niter);
gof = cell(1,Niter);
for it = 1:Niter
    
    % Fit function to fluorescence distribution
    stpt = normrnd(0,1,1,2);
    [fitout{it}, gof{it}] = fit(x(:), y(:), fitfun, 'StartPoint', stpt);
end

% Get parameters from the best iteration
[~,besti] = min(cellfun(@(x) x.sse, gof));
params = coeffvalues(fitout{besti});

% Also return fit output
fitout = fitout{besti};
gof = gof{besti};

end
