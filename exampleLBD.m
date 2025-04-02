%% Load stuff

GPfolder = '/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM';
addpath('/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM')
load(fullfile(GPfolder, 'fullFit.mat'))

%% Part 1: SDoF target design

EALtarget = 0.7;
designSDoF = 1; % update after using the slider

clear options
structureType = 'Wall'; % 'Frame'; % 
options.(structureType).Nstoreys = 4;

lossAssessmentMode = 'global+SLF';
switch lossAssessmentMode
    case 'global' % global level DLRs - structural + nonstructural
        options.FragVuln.damageToLoss = [0 7 15 50 100]/100;

    case 'storey' % storey level DLRs - structural + nonstructural
        options.FragVuln.damageToLoss = [0 1.75 3.75 12.5 25;
                                         0 1.75 3.75 12.5 25;
                                         0 1.75 3.75 12.5 25;
                                         0 1.75 3.75 12.5 25]/100;

    case 'global+SLF'
        options.FragVuln.damageToLoss = [0 7 15 50 100]/100;
        options.NSC.explicitNSC = true;

        % storey loss functions - nonstructural drift [%] sensitive
        SLFshape = 'Papadopoulos';
        NCSDcostRatio = 4/100;
        NSCDpars = [1.66 2.88 1.66 2.89 606.70];
        for n = options.(structureType).Nstoreys : -1 : 1
            options.NSC.Dsens(n).costRatio = NCSDcostRatio;
            options.NSC.Dsens(n).SLFshape = SLFshape;
            options.NSC.Dsens(n).SLFpars = NSCDpars;
        end

        % storey loss functions - nonstructural acceleration sensitive
        NCSAcostRatio = 6/100;
        NSCApars = [2.47 0.79 2.47 0.79 340.54];
        for n = options.(structureType).Nstoreys : -1 : 1
            options.NSC.Asens(n).costRatio = NCSAcostRatio;
            options.NSC.Asens(n).SLFshape = SLFshape;
            options.NSC.Asens(n).SLFpars = NSCApars;
        end
end

options.SDoF.NseedsFy = 60;
options.SDoF.NseedsMu = 60;
options.SDoF.fyBounds = [0.10 0.3];
options.SDoF.minCDR = [0.7 1 1 1]; 

design = lossBasedDesign(structureType, options);
design = DBDgeneral(design);
design = setSeedSDoFproperties(design);
design = getSeedEAL(design, fullFit);
design = getSeedMAFEds(design);
design = getDesignCandidates(design, EALtarget);
design = selectDesignSDoF(design, designSDoF);
plotEAL(design);
plotMAFEds(design);

design = getMemberDetailing(design);
plotMemberDetailing(design)

%% Part 2: design members

% design the frame members to match the requirements of DDBD
