GPfolder = '/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM';
addpath('/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM')
load(fullfile(GPfolder, 'fullFit.mat'))

EALtarget = 0.5;
designSDoF = 147; % update after using the slider

design = lossBasedDesign;
design = getSeedEAL(design, fullFit);
design = getDesignCandidates(design, EALtarget);
design = selectDesignSDoF(design, designSDoF);
plotEAL(design);

design = getFrameMemberDetailing(design);
plotFrameMemberDetailing(design)
