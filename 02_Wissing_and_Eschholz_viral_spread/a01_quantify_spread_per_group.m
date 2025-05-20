% This script estimates the viral spread upon injection of different serotypes and injection volumes as reported in Wissing and Eschholz, PLOS Biology, 2025.  
% The output of each experimental group is saved as a .mat file to be summarized later on. 
%(C)opyright Alexander Dieter, 2025


%% define path and filenames of files to analyze
% filenames = {   'drive:\foldername\examplefile_1';
%                 'drive:\foldername\examplefile_2';
%                 'drive:\foldername\examplefile_3'};


pixelsize = 0.908;      % X/Y-dimension of a single pixel (in µm)
resize_factor = 0.25;   % factor to down-sample image
filterwidth = 50;       % sigma of gaussian filter (in µm)
lvls =  0.1:0.1:0.9;    % define levels of relative fluorescence at which the spread should be compared


viral_spread = [];

for fileIdx = 1:size(filenames, 1)                  % loop across files
    cur_filename = filenames{fileIdx, :};           % get current file

    for hemisphere = 1:2                            % loop across hemispheresto define extension of file name (left/right)

        if hemisphere == 1

            curImage_GFP = imread([cur_filename '_green_left.jpg']);
            curImage_GFP = sum(im2double(curImage_GFP), 3);

        elseif hemisphere == 2
            curImage_GFP = imread([cur_filename '_green_right.jpg']);
            curImage_GFP = sum(im2double(curImage_GFP), 3);

        end

        curImage_GFP_ds = imresize(curImage_GFP,resize_factor);             % re-size image

        filter_sigma = round(filterwidth./(pixelsize./resize_factor));      % calculate width of filter for smoothing
        curImage_GFP_filt = imgaussfilt(curImage_GFP_ds, filter_sigma);     % apply gaussian smoothing to re-sized image
        curImage_GFP_norm = rescale(curImage_GFP_filt, 0, 1);               % normalize brightness


    
        curImage_lvlcontours = zeros(size(curImage_GFP_norm));
        cur_spread_pixels = []; cur_spread_mm2 = [];

        for lvlIdx = 1:length(lvls)    % loop across brightness thresholds as defined above

            curImage_lvlcontours(curImage_GFP_norm > lvls(lvlIdx)) = curImage_lvlcontours(curImage_GFP_norm > lvls(lvlIdx))+10;     % create image based on threshold contours of viral spreads
            cur_spread_pixels(lvlIdx) = sum(sum(curImage_GFP_norm > lvls(lvlIdx)));                                                 % sum up pixels above threshold in order to estimate the viral spread

        end


        %% optional: plot original, downsampled image, filtered and normalized image, and outlines of viral spread at different thresholds.
        %% This might help to find the right parameters for downsampling, smoothing, etc...
%         f0 = figure('color', [1 1 1]);
% 
%         ax1 = subplot(1, 3, 1);
%         imagesc(curImage_GFP_ds);
%         colorbar;  axis image; title('original, downsampled');
%         colormap(ax1,[linspace(0, 0, 500)' linspace(0, 1, 500)' linspace(0, 0, 500)']); 
%         caxis([prctile(curImage_GFP_ds(:), 0.1) prctile(curImage_GFP_ds(:), 99.9)]);
% 
%         ax2 = subplot(1, 3, 2);
%         imagesc(curImage_GFP_norm);
%         colorbar;  axis image; title('filtered and normalized');
%         colormap(ax2,[linspace(0, 0, 500)' linspace(0, 1, 500)' linspace(0, 0, 500)']); 
%     
%         ax3 = subplot(1, 3, 3);
%         imagesc(curImage_lvlcontours);
%         colorbar; axis image; title('topography');
%         colormap(ax3,[linspace(0, 0, 500)' linspace(0, 1, 500)' linspace(0, 0, 500)']); 
% 
%         linkaxes([ax1, ax2, ax3])




        cur_spread_mm2 = (cur_spread_pixels.*(pixelsize./resize_factor).^2)./10^6;      % convert viral spread from pixels to mm²
        viral_spread = vertcat(viral_spread, cur_spread_mm2);                           % collect spread of different samples at different thresholds in one variable

    end

end


savename    =   'test';             % define under which file name data should be stored for subsequent plotting
% save(savename, 'viral_spread')    % uncomment to acutally save data




