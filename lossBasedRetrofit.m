classdef lossBasedRetrofit < lossBasedDesign
    %lossBasedRetrofit
    
    properties
        retrofit
    end


    methods

        function self = lossBasedRetrofit( ...
                structureType, opt, optRetrofit)

            self@lossBasedDesign(structureType, opt);
            self.retrofit = optRetrofit;
            warning('Code some checks for the extra intput')
        end


        function self = DBDgeneral(self)
            
            self = DBDgeneral@lossBasedDesign(self, self.retrofit.Frame.pushBS(2,1));
            
        end


        function self = setSeedSDoFproperties(self)

            % creates the ROWS in the matrices
            fyDummy = linspace(self.parameters.SDoF.fyBounds(1), ...
                self.parameters.SDoF.fyBounds(2), ...
                self.parameters.SDoF.NseedsFy)';
            % creates the COLS in the matrices
            muDummy = linspace(self.parameters.SDoF.muBounds(1), ...
                self.parameters.SDoF.muBounds(2), ...
                self.parameters.SDoF.NseedsMu)';

            [self.FY, self.MU] = meshgrid(fyDummy, muDummy);
            
            % assign MS yield displacement to every seed
            DY = self.retrofit.Frame.deltaYieldMS * ...
                ones(size(self.FY));
            % assign BS yield displacement to the seeds with fy>fy(BS)
            DY(self.FY >= self.retrofit.Frame.pushBS(2,2)) = ...
                self.retrofit.Frame.pushBS(2,1);

            self.T = 2*pi * (DY ./ (self.FY*9.81)).^0.5;

            self.fy = reshape(self.FY, numel(self.FY), 1);
            self.mu = reshape(self.MU, numel(self.MU), 1);
            self.t = reshape(self.T, numel(self.T), 1);

            self.hyst = cell(numel(self.fy),1);
            
            % assign MS hysteresis to the seeds with fy<fy(BS)
            self.hyst(self.fy < self.retrofit.Frame.pushBS(2,2)) = {'MS1'};
            % assign BS hysteresis to the seeds with fy>fy(BS)
            self.hyst(self.fy >= self.retrofit.Frame.pushBS(2,2)) = {'MTf'};

            self.hard = self.parameters.SDoF.hardening * ones(numel(self.fy),1);

            self.fu = self.fy .* (1 + (self.mu-1).*self.hard);
            self.FU = reshape(self.fu, numel(muDummy), numel(fyDummy));
            
            warning('Fix the checkExtrapolation run')
            % checkExtrapolations@lossBasedDesign(self);
            
            % identify indices of the as built structure
            [~, self.retrofit.indicesAsbuilt(1)] =  ...
                min(abs(self.retrofit.Frame.asbuiltPush(3,1) / ...
                self.retrofit.Frame.asbuiltPush(2,1) - self.MU(:,1)));
            [~, self.retrofit.indicesAsbuilt(2)] =  ...
                min(abs(self.retrofit.Frame.asbuiltPush(2,2) - ...
                self.FY(1,:)));

        end


        function self = selectDesignSDoF(self, selectedSDoF, comparativePush)
            
            if nargin < 2; selectedSDoF = 1; end
            if nargin < 3; comparativePush = 0; end

            self = selectDesignSDoF@lossBasedDesign( ...
                self, selectedSDoF, comparativePush);
            plot( get(gca, 'XLim'), ...
                self.retrofit.General.fyFoundation*[1 1], ...
                '--k', 'LineWidth', 2, 'DisplayName', 'Found. Capacity')

            plot(self.retrofit.Frame.asbuiltPush(:,1), ...
                self.retrofit.Frame.asbuiltPush(:,2), ...
                '-', 'Color', [0.6350, 0.0780, 0.1840], ...
                'LineWidth', 2, 'DisplayName', 'As built')
            plot(self.retrofit.Frame.pushBS(:,1), ...
                self.retrofit.Frame.pushBS(:,2) , ...
                '-', 'Color', [0, 0.4470, 0.7410], ...
                'LineWidth', 2, 'DisplayName', 'Beam Sway')

        end


        function plotEAL(self)

            plotEAL@lossBasedDesign(self)
            
            indX = self.retrofit.indicesAsbuilt(1);
            indY = self.retrofit.indicesAsbuilt(2);
            scatter3(self.FY(indX,indY), ...
                self.MU(indX,indY), self.EAL(indX,indY), ...
                1000, 's', 'MarkerFaceColor', [0.64,0.08,0.18], ...
                'DisplayName', 'As built');

            ylim = get(gca, 'Zlim');
            fill3(self.retrofit.General.fyFoundation * ones(4,1), ...
                [get(gca, 'Ylim') fliplr(get(gca, 'Ylim'))]', ...
                [ylim(1) ylim(1) ylim(2) ylim(2)]', ...
                'k', 'FaceAlpha', 0.1, 'DisplayName', 'Found. Capacity')

            fill3(self.retrofit.Frame.pushBS(2,2) *ones(4,1), ...
                [get(gca, 'Ylim') fliplr(get(gca, 'Ylim'))]', ...
                [ylim(1) ylim(1) ylim(2) ylim(2)]', ...
                [0, 0.4470, 0.7410], ...
                'FaceAlpha', 0.1, 'DisplayName', 'Beam Sway')

        end
        

        function plotMemberDetailing(self)
            
            plotMemberDetailing@lossBasedDesign(self)

            plot(self.retrofit.Frame.asbuiltPush(:,1), ...
                self.retrofit.Frame.asbuiltPush(:,2) * ...
                self.parameters.(self.structureType).effMass * 9.81, ...
                '-', 'Color', [0.6350, 0.0780, 0.1840], ...
                'LineWidth', 2, 'DisplayName', 'As built')
            plot(self.retrofit.Frame.pushBS(:,1), ...
                self.retrofit.Frame.pushBS(:,2) * ...
                self.parameters.(self.structureType).effMass * 9.81, ...
                '-', 'Color', [0, 0.4470, 0.7410], ...
                'LineWidth', 2, 'DisplayName', 'Beam Sway')

            legend
        end

        
    end

end