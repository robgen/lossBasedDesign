%% Load stuff

% change this path, clearly :)
GPfolder = '/Users/roberto/Dropbox/2. research/2021 CDR to risk/B_fittingGP/surrogatedPSDM';
addpath(GPfolder)
load(fullfile(GPfolder, 'fullFit.mat'))

%% Input

eqEALtarget = 0.2;
windMAFlimit = 0.02;
designSDoF = 1; % update after using the slider
structureType = 'Frame'; % 'Wall'; % 

opt.SDoF.NseedsFy = 60;
opt.SDoF.NseedsMu = 60;
opt.SDoF.fyBounds = [0.1 0.5];
opt.SDoF.minCDR = [0.7 1 1 1]; 

interstoreyHeight = 4;
opt.Frame.Nparallel = 4;
opt.Frame.Nstoreys = 15;
opt.Frame.massesBuilding = 400*ones(opt.Frame.Nstoreys,1); % [Ton]
opt.Frame.Nbays = 3;

optWind.Building.height = opt.Frame.Nstoreys * interstoreyHeight;
optWind.Building.lengthAcrossWind = 20;
optWind.Building.lengthAlongWind = 20;

optWind.Hazard.timeForAvgWindSpeed = 600; % [s]
optWind.Hazard.densityAir = 1.25; % [kg/m^3]
optWind.Hazard.rugosityTerrain= 1; % [m]
optWind.Hazard.heightDisplacement = 0; % [m] TODO: set name to something sensible
optWind.Hazard.coeffWindward = 0.8; % [-] TODO: check windward is sopravento. set name to something sensible
optWind.Hazard.coeffLeeward = 0.5; % [-] TODO: check leeward is sottovento. set name to something sensible 
optWind.Hazard.typeTerrain = 4; % 1 - 5; Coastal; Open terrain; Sparsely Built-up Suburbs; Towns, Densely Built-up Suburbs; Center of Large Cities; 
optWind.Hazard.accReference = 4; % [cm/s^2] - 6: Offices; 4: Residential
optWind.Hazard.speedWindSeeds = 5:5:50; % [m/s]
optWind.Hazard.indexSpeedWind = 6; % index in the vector "rangeSpeedWind"

optWind.ModellingSAC.betaC = 0.143;
optWind.ModellingSAC.betaD = 0.45;
optWind.ModellingSAC.stDevWeibull = 0.5;

bulkDensity = opt.Frame.massesBuilding(1)*1000 / ...
    (interstoreyHeight * optWind.Building.lengthAcrossWind * optWind.Building.lengthAlongWind);
optWind.SDoF.NseedsBulkDensity = 1;
optWind.SDoF.bulkDensityBounds = [1 1] * bulkDensity;

%% Loss based design

design = LBDearthquakeWind(structureType, opt, optWind);
design = DBDgeneral(design);
design = setSeedSDoFproperties(design);
design = getSeedEAL(design, fullFit);
design = getSeedMAFEds(design);
design = getDesignCandidates(design, eqEALtarget, windMAFlimit);
design = selectDesignSDoF(design, designSDoF);
plotEAL(design);
plotMAFEpsWind(design)
