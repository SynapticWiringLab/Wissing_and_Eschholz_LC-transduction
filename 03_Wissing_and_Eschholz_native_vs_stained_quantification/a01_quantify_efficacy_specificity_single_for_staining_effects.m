% This script compares co-expression of two fluorophores based on the outlines of fluorophore-expressing cells as segmented with CellPose. 
% The script creates a list of files to analyze (based on the content of a defined folder in which data is stored), 
% sequentially loads the data (.tiff-files of microscopy images and .txt
% files defining cell masks), and calculates the number of co-expressing cells based on a pre-defined threshold (overlap).
%(C)opyright Alexander Dieter, 2025


%% clean workspace
clear all; close all; clc;

%% define data
% define the path to the folder in which data is stored. Filenames must end with "_green" for the green channel and "_red" for the red channel. 
% all data files (i.e. tiff-files for both color channels as well as .txt-files from CellPose) should be in this folder, and will be sorted into a list for analysis in the first step of this script. 

path = 'Drive:\datafolder'; 


%% variables to adjust
overlap = 0.5;  % overlap of ROIs to define co-expression



%% create list of files in source folder (i.e. all files that should end up in the analysis)

all_files = dir(fullfile(path,'*.txt'));   
files2analyze = {};                         

for fileIdx = 1:size(all_files)
    curFile = all_files(fileIdx).name;
    
    if contains(curFile,'green')
        dummy = strsplit(curFile,'_green');
        dummy = dummy{1};
    elseif contains( curFile,'red')
        dummy = strsplit(curFile,'_red');
        dummy = dummy{1};
    end
    
    if ~any(strcmp(files2analyze,dummy))
        files2analyze{length(files2analyze)+1} = dummy;
    end
  
end




%% actual analysis
for fileIdx = 1:length(files2analyze)
    curFile = files2analyze{fileIdx};

    clear co_express_THbased co_express_GFPbased

    % load in data
    im_green    = Tiff([path '\' curFile '_green.tif']); im_green = read(im_green);     % load image (green channel)
    im_red      = Tiff([path '\' curFile '_red.tif']); im_red = read(im_red);           % load image (red channel)

    cp_green    = readmatrix([path '\' curFile '_green_cp_outlines.txt']);              % load outlines of detected cells (green channel)
    cp_red      = readmatrix([path '\' curFile '_red_cp_outlines.txt']); 	            % load outlines of detected cells (red channel)

    mastermask_green = f02_create_mastermask(cp_green,  size(im_green));    % create "master mask", i.e. one image containing all segmented cells in a channel (green channel)
    mastermask_red = f02_create_mastermask(cp_red,  size(im_red));          % create "master mask", i.e. one image containing all segmented cells in a channel (red channel)


    [AvBrightRed_THbased, AvBrightGreen_THbased, co_express_THbased] = f03_coexpress(cp_red, mastermask_green, overlap, im_red, im_green);     % calculate co-expression based on red channel 
    [AvBrightGreen_GFPbased, AvBrightRed_GFPbased, co_express_GFPbased] = f03_coexpress(cp_green, mastermask_red, overlap, im_green, im_red);   % calculate co-expression based on green channel 


    dummy = strsplit(curFile, '_'); ID = dummy{1};  % get ID of animal

    % collect data in a commom variable
    d.ID(fileIdx)                      = str2num(ID);                                           % store animal ID for averaging per animal later on 
    d.n_red(fileIdx)                   = length(co_express_THbased);                            % number of detected cells in red channel
    d.n_green(fileIdx)                 = length(co_express_GFPbased);                           % number of detected cells in green channel
    d.co_express{fileIdx}              = co_express_THbased;                                    % TH-based co-expression
    d.co_express_2{fileIdx}            = co_express_GFPbased;                                   % GFP-based co-expression
    d.efficiency(fileIdx)              = sum(co_express_THbased)./length(co_express_THbased);   % efficiency (i.e., TH-based co-expression divided by number of detected TH-positive cells)
    d.specificity(fileIdx)             = sum(co_express_GFPbased)./length(co_express_GFPbased); % specificity (i.e., GFP-based co-expression divided by number of detected GFP-positive cells)
    
     d.AvBrightRed_THbased{fileIdx}          = AvBrightRed_THbased;                             % average brightness of red channel in TH-based cell masks            
    d.AvBrightGreen_THbased{fileIdx}        = AvBrightGreen_THbased;                            % average brightness of green channel in TH-based cell masks            
    d.AvBrightRed_GFPbased{fileIdx}         = AvBrightRed_GFPbased;                             % average brightness of red channel in GFP-based cell masks            
    d.AvBrightGreen_GFPbased{fileIdx}       = AvBrightGreen_GFPbased;                           % average brightness of green channel in GFP-based cell masks            

    disp([num2str(fileIdx) ' out of ' num2str(length(files2analyze)) ' images analyzed']) % provide feedback on progress of protocol


end



%% save data for subsequent quantification
savename = 'test';      % define the filename in which your data should be saved
% save(savename, 'd');   % uncomment  to save data 
