classdef LBDearthquakeWind < lossBasedDesign
    %LBDearthquakeWind
    
    properties
        windParams

        freqFirstMode
        bulkDensity
        MAFdisplacement
        MAFacceleration
        isCandidateEQ_W
    end


    methods

        function self = LBDearthquakeWind( ...
                structureType, opt, optWind)

            self@lossBasedDesign(structureType, opt);
            self.windParams = optWind;
            warning('Code some checks for the extra input')
        end


        function self = getDesignCandidates(self, eqEALtarget, windMAFlimit)
            % earthquake-based candidates
            self = getDesignCandidates@lossBasedDesign(self, eqEALtarget);
            
            % earthquake and wind-based candidates
            self = setWindSeedProperties(self);
            self = getWindSeedMAFEps(self, windMAFlimit);
            self = getDesignCandidatesWind(self, windMAFlimit);
        end


        function self = setWindSeedProperties(self)
            
            numEarthquakeCandidates = sum(sum(self.isCandidate));
            n1EarthquakeCandidates = 1 ./ self.T(self.isCandidate);
            seedBulkDensity = ceil(linspace( ...
                self.windParams.SDoF.bulkDensityBounds(1), ...
                self.windParams.SDoF.bulkDensityBounds(2), ...
                self.windParams.SDoF.NseedsBulkDensity)); % TODO: a sprintf only works if this vector has integers in it. Fix it there

            self.freqFirstMode = repmat( n1EarthquakeCandidates, ...
                1, self.windParams.SDoF.NseedsBulkDensity);
            self.bulkDensity = repmat( ...
                seedBulkDensity, numEarthquakeCandidates, 1);

        end


        function self = getWindSeedMAFEps(self, windMAFlimit)

            sys.h = self.windParams.Building.height;
            sys.b = self.windParams.Building.lengthAcrossWind;
            sys.d = self.windParams.Building.lengthAlongWind;

            wind.T = self.windParams.Hazard.timeForAvgWindSpeed;
            wind.rho_air = self.windParams.Hazard.densityAir;
            wind.z_0 = self.windParams.Hazard.rugosityTerrain;
            wind.z_d = self.windParams.Hazard.heightDisplacement;
            wind.C_w = self.windParams.Hazard.coeffWindward;
            wind.C_l = self.windParams.Hazard.coeffLeeward;
            wind.Type_of_terrain = self.windParams.Hazard.typeTerrain;
            wind.a_0 = self.windParams.Hazard.accReference;
            rangeVento = self.windParams.Hazard.speedWindSeeds;

            beta_c = self.windParams.ModellingSAC.betaC;
            beta_d = self.windParams.ModellingSAC.betaD;
            dev_std_norm_WEI = self.windParams.ModellingSAC.stDevWeibull;
            
            MAF_limite_acc = windMAFlimit;
            MAF_limite_spost = windMAFlimit;
            
            rho_bulk = self.bulkDensity(1,:);
            n1 = self.freqFirstMode(:,1)';

            [MAF_disp, MAF_acc] = SACFEMA_wind_fun( ...
                rho_bulk, n1, rangeVento, ...
                MAF_limite_acc, MAF_limite_spost, ...
                sys, wind, beta_d, beta_c, dev_std_norm_WEI);

            referenceWindSpeed = self.windParams.Hazard.speedWindSeeds( ...
                self.windParams.Hazard.indexSpeedWind);
            self.MAFacceleration = MAF_acc.( ...
                sprintf('Vel_%02d', referenceWindSpeed))';
            self.MAFdisplacement = MAF_disp.( ...
                sprintf('Vel_%02d', referenceWindSpeed))';
            
        end


        function self = getDesignCandidatesWind(self, windMAFlimit)

            self.isCandidateEQ_W = self.MAFdisplacement<windMAFlimit & ...
                self.MAFacceleration<windMAFlimit;

        end


        function self = selectDesignSDoF(self, selectedSDoF, comparativePush)
            
            if nargin < 2; selectedSDoF = 1; end
            if nargin < 3; comparativePush = 0; end

            self = selectDesignSDoF@lossBasedDesign( ...
                self, selectedSDoF, comparativePush);
            
        end


        function plotEAL(self)

            plotEAL@lossBasedDesign(self)

        end
        

        function plotMAFEpsWind(self)
            
            colPS1 = [0.000 0.447 0.741];
            colPS2 = [0.850 0.325 0.098];

            figure;
            scatter3(self.freqFirstMode(self.isCandidateEQ_W), ...
                self.bulkDensity(self.isCandidateEQ_W), ...
                self.MAFdisplacement(self.isCandidateEQ_W), ...
                50, colPS1, 'filled', 'DisplayName', 'PS1 (displacement)')
            hold on;
            scatter3(self.freqFirstMode(self.isCandidateEQ_W), ...
                self.bulkDensity(self.isCandidateEQ_W), ...
                self.MAFacceleration(self.isCandidateEQ_W), ...
                30, colPS2, 'filled', 'DisplayName', 'PS2 (acceleration)')
            
            scatter3(self.freqFirstMode(~self.isCandidateEQ_W), ...
                self.bulkDensity(~self.isCandidateEQ_W), ...
                self.MAFdisplacement(~self.isCandidateEQ_W), ...
                50, colPS1, 'DisplayName', 'PS1 - Non comp.')
            scatter3(self.freqFirstMode(~self.isCandidateEQ_W), ...
                self.bulkDensity(~self.isCandidateEQ_W), ...
                self.MAFacceleration(~self.isCandidateEQ_W), ...
                30, colPS2, 'DisplayName', 'PS2 - Non comp.')

            set(gca, 'FontSize', 18)
            xlabel('First mode frequency [Hz]')
            ylabel('Building bulk density [kg/m3]')
            zlabel('MAFE performance state [1/y]')
            legend1 = legend;
            set(legend1, 'Position',[0.583 0.716 0.338 0.199]);
            view([47.421602641144 12.6164650749273]);
            
        end
        

    end

end