%% Load stuff

GPfolder = '/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM';
addpath('/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM')
load(fullfile(GPfolder, 'fullFit.mat'))

%% Part 1: SDoF target design

EALtarget = 0.5;
designSDoF = 117; % update after using the slider

design = lossBasedDesign;
design = getSeedEAL(design, fullFit);
design = getDesignCandidates(design, EALtarget);
design = selectDesignSDoF(design, designSDoF);
plotEAL(design);

design = getFrameMemberDetailing(design);
plotFrameMemberDetailing(design)

%% Part 2: design members

% design the frame members to match the requirements of DDBD
% create a slamasolver file
% create the ruaumoko model from SLaMAsolver

%% Part 3: 

