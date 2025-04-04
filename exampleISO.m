%%  Example of DLBD of isolated structure using LRB or FPS
% clear; clc
%% LOAD SURROGATE PSDMs
GPfolder = 'C:\Users\dsuas\My_Drive\Diego Suarez DLBD\3. Codes\surrogatedPSDM_IsolatedStructures_BU'; %add folder where GP regression is located. Refer to github repo surrogatePSDMsISO
addpath('C:\Users\dsuas\My_Drive\Diego Suarez DLBD\3. Codes\surrogatedPSDM_IsolatedStructures_BU')
addpath('C:\Users\dsuas\My_Drive\Diego Suarez DLBD\3. Codes\lossBasedDesign')
addpath('C:\Users\dsuas\My_Drive\Diego Suarez DLBD\3. Codes\Roberto_Codes\robSeismicAnalyses')
load(fullfile(GPfolder, 'LRBfullFitDisp.mat'))
load(fullfile(GPfolder, 'LRBfullFitAcc.mat'))

%%   GENERAL PARAMETERS AND INPUT:  
% clearvars -except fullFit fullFitAcc GPfolder designC1

% CONTROLS
options.Controls.DET=1; %1 for defining the solution space in terms of 
                        %the isolation system characteristics
options.Controls.ColLosses=0; %1 for considering collapse losses

% STRUCTURE TYPE
isolationtype = 'LRB'; 
sstype='Wall';

% design REQUIREMENTS
options.General.EALtarget = 0.3;
options.General.toleranceEAL = 0.005;
options.General.MAFECollapseLim = 1.2E-5; % 3% given MCE
options.General.MAFESSYieldLim = 5E-5; 
options.General.lossType = 'repairTime';

options.LossParms.numWorkersPerFloor=1; 
options.LossParms.workingHoursRatio=1;
options.LossParms.workDaysRatio= 1;
options.LossParms.repairScheme='series'; %Parallel or series
% options.LossParms.indLossesPerDay=100000;

% LEAD RUBBER BEARINGS SYSTEM PARAMETERSso
% Specify wanted seed structure's iso system's bounds for 
% Height, ARatio and Equivalent Axial Stress
options.Iso.costRatio = 1;
options.Iso.DS = [0.4 1]; % Fraction of ultimate ductility [-]
options.Iso.DLR = [0 3 584]; % Last DLR corresponds to collapse cost for entire building
options.Iso.Mu= 13.33; % Ultimate displacement ductility
options.Iso.HeightBounds=[0.15 0.4]; % [m]
options.Iso.ARatioBounds=[8 20]; % Rubber to lead areas ratio [-]
options.Iso.SigmaBounds=[6 20]; % Total seismic weight of the structure divided by total bearing area of isolators [MPa]
options.Iso.NseedsHeight = 10;
options.Iso.NseedsARatio = 10;
options.Iso.NseedsSigma = 10;
options.Iso.Glead=130; % Shear modulus of lead [MPa]
options.Iso.Grubber=1; % Shear modulus of rubber [MPa]
options.Iso.epsYlead=0.075; % Yield strain of lead 
options.Iso.fylead=10; % Yield stress of lead [MPa]

% SUPERSTRUCTURE PARAMETERS (SS):
options.SS.costRatio = 1;
options.SS.DS = [0.5 1];
options.SS.DLR = [0 6 20];
options.SS.FyBounds = [0.18 0.18];
options.SS.NseedsFy = 1;
options.SS.FloorsWeight=[6630 7500 7500  7500]; %Include the weight of the isolation base on the first entry [kN] 
options.SS.InterStoreyHeights=[0.5 3 3 3];%Include the height of the CM of the isolation level [m]
options.SS.numStoreys=length(options.SS.FloorsWeight)-1;
options.SS.partFactor=0.7514;
options.SS.meffFactor=0.7514;
options.SS.DYeff=0.025; %[m] Optional
options.SS.DYatTop=0.0332; %[m] Optional 
options.SS.t1=2*pi.*(options.SS.meffFactor./(options.SS.FyBounds/options.SS.DYeff*9.81)).^0.5; %Fundamental period considering the non-isolated configuration [s] 
options.SS.accProf='Linear'; 




%%
lossAssessmentMode = 'global+SLF';
switch lossAssessmentMode
    
    case 'global+SLF'

        options.NSC.explicitNSC = true;
        
        % storey loss functions - nonstructural drift [%] sensitive
        SLFshape = 'Weibull';
        NSCDpars = [0.97653728 2.18239469 2.44736503];
        NSCDrepairTime = [116 116 116]; % same for all stories
        for n = options.SS.numStoreys : -1 : 1
            options.NSC.Dsens(n).totalRepairTime = NSCDrepairTime(n);
            options.NSC.Dsens(n).SLFshape = SLFshape;
            options.NSC.Dsens(n).SLFpars = NSCDpars;
        end

        % storey loss functions - nonstructural acceleration [g] sensitive
        NSCApars = [  0.97492696 0.68261267 1.79586322;
                    0.98135802 1.57794117 2.13571581;
                    0.98135802 1.57794117 2.13571581;
                    0.97920735 1.81444311 3.57037361];
        NSCArepairTime = [185 460 460 326]; 
        for n = options.SS.numStoreys+1 : -1 : 1
            options.NSC.Asens(n).totalRepairTime = NSCArepairTime(n);
            options.NSC.Asens(n).SLFshape = SLFshape;
            options.NSC.Asens(n).SLFpars = NSCApars(n,:);
        end
end

%%
% FRAGILITY AND VULNERABILITY PARAMETERS
options.FragVuln.fixedBeta=NaN; 
options.FragVuln.betaSDoFtoMDoF=0.4;
options.FragVuln.maxIM=3.5; % [g]
options.FragVuln.samplesIM=1300;

% SITE SPECIFIC SEISMIC HAZARD
options.Hazard.faultRate=0.58;
options.Hazard.periodsHazard=[0,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1,1.5,2];
options.Hazard.extrapMode='powerLaw';
options.Hazard.MAFEhazard=[0.033214624,0.019885045,0.013862944,0.009885926,0.007133499,0.004969227,0.00210721,0.001025866,0.000404054];
options.Hazard.intensityHazard=[0.118500000000000,0.231901076160571,0.284400000000000,0.284400000000000,0.284400000000000,0.284400000000000,0.249647374490702,0.166427989012818,0.124820042771853,0.0832129082285958,0.0597891517804315;0.156000000000000,0.294972075787921,0.363792000000000,0.363792000000000,0.363792000000000,0.363792000000000,0.326382422300220,0.217582434028613,0.163185282252471,0.108789450714589,0.0815918939267367;0.183000000000000,0.341305827038763,0.420458740558145,0.424194000000000,0.424194000000000,0.424194000000000,0.387811782857298,0.258529430931683,0.193893977166246,0.129261174159895,0.0969454922187405;0.213000000000000,0.392400522703245,0.482100784054867,0.490752000000000,0.490752000000000,0.490752000000000,0.455921264959425,0.303930770035728,0.227943687239298,0.151960368720452,0.113969728298491;0.241667382400000,0.438980114035280,0.537636479852921,0.556076646902400,0.556076646902400,0.556076646902400,0.531683535429199,0.354443451456181,0.265829361119539,0.177218032623064,0.132913118895613;0.274027891000000,0.495860969580310,0.606777508870465,0.634374567665000,0.634374567665000,0.634374567665000,0.618336644289972,0.412206673094016,0.309150327439777,0.206097985927785,0.154572902412411;0.347077173600000,0.621962545909488,0.759405232064231,0.820490438390400,0.820490438390400,0.820490438390400,0.820490438390400,0.565248464301760,0.423929063051237,0.282615900631437,0.211961012200778;0.407159360000000,0.727702500307648,0.887974070461473,0.977182464000000,0.977182464000000,0.977182464000000,0.977182464000000,0.695119133694674,0.521330086813096,0.347548986783303,0.260660585294040;0.467092460800000,0.836572020006089,1.02131179960913,1.14811326864640,1.14811326864640,1.14811326864640,1.14811326864640,0.846507595616605,0.634870338466778,0.423241969970881,0.317430187426319];

%% DLBD_IsolatedStructures Execution

design = lossBasedDesignISO(isolationtype, sstype,...
     options);
design = getSSparameters(design,options.SS.DYatTop);
design = getSeedEAL(design, fullFit, fullFitAcc);
design = getSeedMAFEds(design);
design = getDesignCandidates(design);
design = selectDesign(design,[696 1]); % [indISO,indSS]
design = getDesignProperties(design);
plotCandidates(design);
% clearvars -except design
