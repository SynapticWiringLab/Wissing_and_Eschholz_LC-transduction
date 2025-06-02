# Wissing_and_Eschholz_LC-transduction

This repository contains custom MATLAB-code that allows reproduction of the analysis of data presented in the manuscript 
"A comparison of viral strategies and model systems to target norepinephrine neurons in the locus coeruleus reveals high variability in transgene expression patterns" 
by Wissing and Eschholz et al., published in PLOS Biology, 2025 (doi: 10.1371/journal.pbio.3003228). 

The corresponding data is deposited at: 
https://gin.g-node.org/SW_lab/Wissing_and_Eschholz_et_al_PLOS_Biology_2025/ (doi: 10.12751/g-node.msn52q)

The path linking to the data has to be adjusted in the beginning of each script (which is commented in the individual scripts). 
The code was tested on MATLAB version 2022b.


This repository contains the following analyses of the paper: 

(1) Quantification of the efficacy and specificity of transgene expression of immunostained coronal sections of the mouse Locus coeruleus, 
     as well as an analysis of the relative brightness of fluorescence. The script to process data (maximum projections in .tif format along with CellPose segmented cell masks in .txt format) #
     is called "a01_quantify_efficacy_specificity_single". The following scripts are used to plot the data resulting from this script: 
          a) efficacy of transgene expression across model systems (Fig. 1, panel D): a02_plot_Fig1_panelD_efficacy
          b) specificty of transgene expression across model systems (Fig. 1, panel E): a02_plot_Fig1_panelE_specificity
          c) Relative eGFP brightness in TH-eGFP-co-expressing cells vs cells exclusively expressing eGFL (Fig. 1, panel F): a02_plot_Fig1_panelF_GFPlevels
          d) efficacy of transgene expression using CAG-DIO-eGFP vs hSyn-DIO-eGFP in DBH-cre mice (Fig. 3, panel C): a02_plot_Fig3_panelC_efficacy_CAG_vs_hSyn
          e) specificity of transgene expression using CAG-DIO-eGFP vs hSyn-DIO-eGFP in DBH-cre mice (Fig. 3, panel C): a02_plot_Fig3_panelD_specificity_CAG_vs_hSyn

To reproduce our results, one can either run all scripts using the raw input data provided on the g-node repository specified above, or run the scripts a-e for visualization of the data, which is already provided in the github-repository (*50_WissingEschholz2025.mat)


(2) Quantification of the spatial extent of viral spread upon rAAV-based neuronal transduction in the mouse LC. 
    The scripts to quantify the spatial spread per experimental group is called "a01_quantify_spread_per_group". The following scripts are used to plot the data resulting from this script: 
          a) viral spread by serotype, comparing 300 nl of rAAV2/9 vs 300 nl of rAAV2/2 (Fig. 3, panel G): a02_plot_Fig3_panelG_viral_spread_by_serotype
          b) viral spread by injection volume, comparing 50, 100, and 300 nl of rAAV2/9 (Fig. 3, panel J): a02_plot_Fig3_panelJ_viral_spread_by_injection_volume

To reproduce our results, one can either run all scripts using the raw input data provided on the g-node repository specified above, or run the scripts a-e for visualization of the data, which is already provided in the github-repository (*WissingEschholz2025.mat)


(3) The effect of immuno-staining as compared to native fluorescent brain slices on subsequenta nalysis reported in the manuscript (Supplementary Figure S7). 
    The script to process data (maximum projections in .tif format along with CellPose segmented cell masks in .txt format) #
    is called "a01_quantify_efficacy_specificity_single_for_staining_effects" and is very similar to the script in (1).
    Output data is then plotted by using script "a02_LCt_FigureS7_20231016".

To reproduce our results, one can either run all scripts using the raw input data provided on the g-node repository specified above, or run the scripts a-e for visualization of the data, which is already provided in the github-repository (*50_WissingEschholz2025.mat) 

    
    

