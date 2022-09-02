%% Load stuff

% change this path, clearly :)
GPfolder = '/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM';
addpath(GPfolder)
load(fullfile(GPfolder, 'fullFit.mat'))

%% Part 1: SDoF target design

EALtarget = 0.7;
designSDoF = 1; % update after using the slider
structureType = 'Frame'; % 'Wall'; % 

options.SDoF.NseedsFy = 60;
options.SDoF.NseedsMu = 60;
options.SDoF.fyBounds = [0.10 0.3];
options.SDoF.minCDR = [0.7 1 1 1]; 

design = lossBasedRetrofit(structureType, options);
design = getSeedEAL(design, fullFit);
design = getSeedMAFEds(design);
design = getDesignCandidates(design, EALtarget);
design = selectDesignSDoF(design, designSDoF);
plotEAL(design);

design = getMemberDetailing(design);
plotMemberDetailing(design)

%% Part 2: design members

% design the frame members to match the requirements of DDBD
