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
            self = setSeedPropertiesWind(self);
            self = getSeedMAFEpsWind(self, windMAFlimit);
            self = getDesignCandidatesWind(self, windMAFlimit);
        end


        function self = setSeedPropertiesWind(self)
            
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


        function self = getSeedMAFEpsWind(self, windMAFlimit)

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

            [MAF_disp, MAF_acc] = self.getMAFEwindPSs( ...
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

    methods(Static)
        
        function [lambda_PS_spost,lambda_PS_acc] = getMAFEwindPSs( ...
                rho_bulk, n1, rangeVento, MAF_limite_acc, MAF_limite_spost, ...
                sys, wind, beta_d, beta_c, dev_std_norm_WEI)

            Range.T_ritorno = [1 50];       % Periodo di ritorno per la velocità del vento [anni]
            % 1: Accelerazione
            % 50: Spostamento
            Range.Vento = rangeVento;           % Velocità di base per il vento [m/s]: Ho bisogno di più velocità per costruire delle curve da interpolare
            % 10: Bassa ventosità
            % 30: Alta ventosità
            Range.rho_bulk = rho_bulk;      % Bulk density [kg/m^3] - Ipotesi che sia costante con l'altezza
            Range.n1 = n1;                  % Frequenza fondamentale along-wind [Hz]

            for t = 1:length(Range.T_ritorno)
                for v = 1:length(Range.Vento)
                    for i = 1:length(Range.rho_bulk)
                        for j = 1:length(Range.n1)


                            %%  DEFINIZIONE DATI
                            %Questi parametri sono ricavati analiticamente da i dati inseriti in
                            %precedenza quindi   ----->   NON VANNO CAMBIATI

                            if Range.T_ritorno(t) == 50         % SLU spostamento
                                sys.xi_1 = 0.012;   % Coefficiente di smorzamento modo along-wind [-]
                            elseif Range.T_ritorno(t) == 1      % SLE accelerazione
                                sys.xi_1 = 0.008;   % Coefficiente di smorzamento modo along-wind [-]
                            end
                            sys.n1 = Range.n1(j);   % Frequenza fondamentale along-wind [Hz]
                            sys.rho_bulk = Range.rho_bulk(i);   % Bulk density [kg/m^3]

                            wind.T_ritorno = Range.T_ritorno(t);
                            wind.U_base = Range.Vento(v);     % Velocità di base a z=10m [m/s]

                            if wind.T_ritorno == 1
                                wind.cr = 0.75;
                            elseif wind.T_ritorno == 50
                                wind.cr = 1;
                            else
                                error('Return period can only be 1 or 50')
                            end

                            wind.U_ref = wind.U_base * wind.cr;     % Velocità di riferimento a z=10m [m/s]
                            wind.U_ref = wind.U_base;     % Velocità di riferimento a z=10m [m/s]

                            switch wind.Type_of_terrain
                                case 1
                                    wind.beta = 6.5;
                                    wind.u_star_ratio = 0.85;
                                case 2
                                    wind.beta = 6.0;
                                    wind.u_star_ratio = 1.00;
                                case 3
                                    wind.beta = 5.25;
                                    wind.u_star_ratio = 1.15;
                                case 4
                                    wind.beta = 4.85;
                                    wind.u_star_ratio = 1.33;
                                case 5
                                    wind.beta = 4.00;
                                    wind.u_star_ratio = 1.45;
                                otherwise
                                    error('Invalid value for T. Must be between 1 and 5.');
                            end

                            wind.u_star_ref = wind.U_ref / (2.5 * log((10 - wind.z_d) / 0.07));     % Reference Fiction velocity [m/s]
                            wind.u_star = wind.u_star_ref * wind.u_star_ratio;                      % Fiction velocity [m/s]
                            wind.C_d = wind.C_w + wind.C_l;         % Coefficiente di drag [-]

                            %% RISULTATI

                            % Solari 1982 - Along-wind response estimation
                            [Spost_max, Acc_max] = LBDearthquakeWind.get_top_disp_acc_Solari(sys, wind);

                            FieldName = sprintf('TRitorno_%02d_Vento_%02d', Range.T_ritorno(t), Range.Vento(v));
                            if Range.T_ritorno(t) == 50         % SLU spostamento
                                Risultati.(FieldName)(i,j) = Spost_max;
                            elseif Range.T_ritorno(t) == 1
                                Risultati.(FieldName)(i,j) = Acc_max;
                            end

                            if sys.n1 < 1
                                n_limite = wind.a_0 / (sys.n1^0.56);
                            elseif (1 <= sys.n1) && (sys.n1 <=2)
                                n_limite = wind.a_0;
                            else
                                n_limite = 0.5 * wind.a_0 * sys.n1;
                            end

                            Acc_limite(i,j) = n_limite / 100;                       % Per passare da cm/s^2 a m/s^2
                            Spost_limite(i,j) = 80 / 500;


                            FieldName2 = sprintf('TRitorno_%02d_Bulk_%03d_freq_%03d', Range.T_ritorno(t), Range.rho_bulk(i), round(Range.n1(j)*100));
                            Curve_D_IntMeas.(FieldName2)(2,v) = Risultati.(FieldName)(i,j);         %Curve da interpolare che mi servono dopo per SAC-FEMA
                            Curve_D_IntMeas.(FieldName2)(1,v) = wind.U_base;
                        end
                    end
                end
            end


            %% Procedura SAC-FEMA PER VENTO

            for t = 1:length(Range.T_ritorno)
                for v = 1:length(Range.Vento)
                    for i = 1:length(Range.rho_bulk)   %8
                        for j = 1:length(Range.n1)   %40


                            if Range.T_ritorno(t) == 50         % Deviazione standard della domanda
                                beta_d = 0.273;
                            elseif Range.T_ritorno(t) == 1
                                beta_d = 0.449;
                            end


                            [a_i,b_i] = LBDearthquakeWind.Interp_1(Range,Curve_D_IntMeas,i,j,t);

                            % plot(x,Y,'ko',x,(a_i*x.^b_i)/1000,'b-')

                            [Vettore_V10] = LBDearthquakeWind.Vettore_vel_10(Range,v,a_i,b_i,beta_d,beta_c);

                            % Interpolazione per trovare i tre valori di H(im)

                            [H_values, k, lambda] = LBDearthquakeWind.interpolating_law( ...
                                Range.Vento(v), dev_std_norm_WEI, Vettore_V10);
                            % x = 1:1:50;
                            % plot(x, 1 - wblcdf(x, lambda, k))

                            [k_sol] = LBDearthquakeWind.solve_k(H_values, Vettore_V10);

                            [MAF, FieldName3] = LBDearthquakeWind.Calcolo_MAF(Range, ...
                                Spost_limite, Acc_limite, k_sol, beta_d,...
                                b_i, beta_c, i, j, t, v);

                            if Range.T_ritorno(t) == 50         % SLU spostamento
                                lambda_PS_spost.(FieldName3)(i,j) = MAF;
                            elseif Range.T_ritorno(t) == 1
                                lambda_PS_acc.(FieldName3)(i,j) = MAF;
                            end

                            %nella matrice lambda_PS_ è associata un valore della MAF per ogni coppia
                            %di punti frequenza, bulk density a seconda che stiamo parlando di
                            %accelerazione o spostamento

                            MAF_limite_acc(i,j) = 0.02;
                            MAF_limite_spost(i,j) = 0.02;

                        end
                    end
                end
            end

        end


        function [x_peak,xddot_peak] = get_top_disp_acc_Solari(sys,wind)

            % method by Solari TODO: specify which paper
            Calcolo.Q = 2 * (1 - wind.z_d^2 / sys.h^2) * log((sys.h - wind.z_d) / wind.z_0) - 1;
            Calcolo.J = 0.78 * Calcolo.Q^2;
            Calcolo.BB = 6.71 * Calcolo.Q^2 / (1 + 0.26 * sys.b / sys.h);
            Calcolo.f1 = sys.n1 * sys.h / wind.u_star;
            Calcolo.C_x = 0;
            Calcolo.x_1 = 12.32 * Calcolo.f1 * sys.d / (Calcolo.Q * sys.h);
            Calcolo.N_f1 = 1 / Calcolo.x_1 - 1 / (2 * Calcolo.x_1^2) * (1 - exp(-2 * Calcolo.x_1));
            Calcolo.C_Df_f1 = sqrt(wind.C_w^2 + 2 * wind.C_w * wind.C_l * Calcolo.N_f1 + wind.C_l^2);
            Calcolo.x_2 = 3.55 * Calcolo.f1 / Calcolo.Q;
            Calcolo.M_z = sys.b * sys.d * sys.rho_bulk;
            Calcolo.C_x2 = 1 / Calcolo.x_2 - 1 / (2 * Calcolo.x_2^2) * (1 - exp(-2 * Calcolo.x_2));
            Calcolo.RR = 0.59 * Calcolo.Q^2 / sys.xi_1 * (Calcolo.Q / Calcolo.f1)^(2/3) * ...
                Calcolo.C_Df_f1^2 / wind.C_d^2 * Calcolo.C_x2 * ...
                1 / (1 + 3.95 * (Calcolo.f1 /Calcolo.Q) * (sys.b / sys.h));
            Calcolo.M_1 = Calcolo.M_z * sys.h / 3;
            Calcolo.q_star = 0.5 * wind.rho_air * wind.u_star^2;
            Calcolo.x_mean = wind.C_d * sys.b * sys.h * Calcolo.q_star * Calcolo.J / ...
                (Calcolo.M_1 * (2 * pi * sys.n1)^2);
            Calcolo.sigma_x = wind.C_d * sys.b * sys.h * Calcolo.q_star / ...
                (Calcolo.M_1 * (2 * pi * sys.n1)^2) * sqrt(wind.beta * Calcolo.BB / 6 + Calcolo.RR);
            Calcolo.v_x = sys.n1 * sqrt(Calcolo.RR / (wind.beta * Calcolo.BB / 6 + Calcolo.RR));
            Calcolo.K_x = (1.175 + 2 * log(Calcolo.v_x * wind.T))^(1/2);
            Calcolo.G = 1 + Calcolo.K_x * Calcolo.sigma_x / Calcolo.x_mean;
            Calcolo.X_max = Calcolo.G * Calcolo.x_mean;
            Calcolo.sigma_xddot = wind.C_d * sys.b * sys.h * Calcolo.q_star / Calcolo.M_1 * sqrt(Calcolo.RR);
            Calcolo.K_xddot = (1.175 + 2 * log(sys.n1 * wind.T))^(1/2);
            Calcolo.Xddot_max = Calcolo.K_xddot * Calcolo.sigma_xddot;

            x_peak = Calcolo.X_max;
            xddot_peak = Calcolo.Xddot_max;

        end


        function [a_i,b_i] = Interp_1(Range,Curve_D_IntMeas,i,j,t)
            % Funzione per interpolazione della curva parametro intensità vs prestazione media


            FieldName2 = sprintf('TRitorno_%02d_Bulk_%03d_freq_%03d', Range.T_ritorno(t), Range.rho_bulk(i), round(Range.n1(j)*100));

            % Interpolazione per a_i e b_i

            x = Curve_D_IntMeas.(FieldName2)(1,4:length(Range.Vento));
            Y = Curve_D_IntMeas.(FieldName2)(2,4:length(Range.Vento));

            log_x = log(x);
            log_Y = log(Y * 1000);
            p = polyfit(log_x, log_Y, 1);
            b_i = p(1);
            a_i = exp(p(2));
            % a_i = 0.00056;


        end
    

        function Vettore_V10 = Vettore_vel_10(Range, v, a_i, b_i, beta_d, beta_c)
            % Funzione per interpolazione della curva parametro intensità vs prestazione media


            Vettore_V10(1,1) = (Range.Vento(v) / a_i) ^ (1/b_i) * ...
                exp(-0.5 * sqrt(beta_d^2 + beta_c^2) / b_i);

            Vettore_V10(2,1) = (Range.Vento(v) / a_i) ^ (1/b_i) * ...
                exp(-1.5 * sqrt(beta_d^2 + beta_c^2) / b_i);

            Vettore_V10(3,1) = (Range.Vento(v) / a_i) ^ (1/b_i) * ...
                exp(-3.0 * sqrt(beta_d^2 + beta_c^2) / b_i);
        end


        function [H_values, k, lambda] = interpolating_law(mu_WEI, dev_std_norm_WEI, Vettore_V10)
            % Utilizzo della distribuzione di weibull per la velocità del vento

            sigma = dev_std_norm_WEI * mu_WEI;

            % Funzione anonima per trovare k e lambda
            f = @(p) [p(2) * gamma(1 + 1/p(1)) - mu_WEI; ...
                p(2)^2 * (gamma(1 + 2/p(1)) - (gamma(1 + 1/p(1)))^2) - sigma^2];

            % Valori iniziali (stima approssimativa)
            p0 = [1, 1];

            % Risoluzione numerica
            options = optimoptions('fsolve', 'Display', 'off');
            [k_lambda, ~] = fsolve(f, p0, options);

            k = k_lambda(1); % Forma
            lambda = k_lambda(2); % Scala

            H_values = 1 - wblcdf(Vettore_V10, lambda, k);

        end
    
    
        function k_values = solve_k(H_values, Vettore_V10)
            % Funzione per trovare i valori di k0, k1 e k2
            %
            % H_values: vettore con i tre valori noti di H
            % s_values: vettore con i tre valori noti di s

            % % % Funzione anonima per il sistema di equazioni
            % % fun = @(k) [
            % %     k(1) * exp(-k(3) * (log(Vettore_V10(1)))^2 - k(2) * log(Vettore_V10(1))) - H_values(1);
            % %     k(1) * exp(-k(3) * (log(Vettore_V10(2)))^2 - k(2) * log(Vettore_V10(2))) - H_values(2);
            % %     k(1) * exp(-k(3) * (log(Vettore_V10(3)))^2 - k(2) * log(Vettore_V10(3))) - H_values(3)
            % % ];
            % %
            % % % Valori iniziali per fsolve
            % % k0_guess = [1e-18; -30; 8];
            % %
            % % % Risoluzione del sistema non lineare
            % % options = optimoptions('fsolve', 'Display', 'iter');
            % % k_values = fsolve(fun, k0_guess, options);

            k1 = (log(H_values(3)/H_values(1)) + log(H_values(2)/H_values(1)) * (log(Vettore_V10(3)))^2 / ((log(Vettore_V10(1)))^2 - (log(Vettore_V10(2)))^2)...
                - log(H_values(2)/H_values(1)) * (log(Vettore_V10(1)))^2 / ((log(Vettore_V10(1)))^2 - (log(Vettore_V10(2)))^2)) / ...
                ((log(Vettore_V10(1)) - log(Vettore_V10(2))) * 1 / ((log(Vettore_V10(1)))^2 - (log(Vettore_V10(2)))^2) * (log(Vettore_V10(3)))^2 - log(Vettore_V10(3)) - ...
                (log(Vettore_V10(1)) - log(Vettore_V10(2))) * 1 / ((log(Vettore_V10(1)))^2 - (log(Vettore_V10(2)))^2) * (log(Vettore_V10(1)))^2 + log(Vettore_V10(1)));

            k2 = (log(H_values(2)/H_values(1)) - k1 * (log(Vettore_V10(1)) - log(Vettore_V10(2)))) * 1 / ((log(Vettore_V10(1)))^2 - (log(Vettore_V10(2)))^2);

            k0 = H_values(1) / exp(-k2 * (log(Vettore_V10(1)))^2 - k1 * log(Vettore_V10(1)));

            k_values(1) = k0;
            k_values(2) = k1;
            k_values(3) = k2;

        end
    
    
        function [MAF,FieldName3] = Calcolo_MAF( ...
                Range, Spost_limite, Acc_limite, k_sol, beta_d,...
                b_i, beta_c, i, j, t, v)
                % Funzione per interpolazione della curva parametro intensità vs prestazione media


            q_i = 1 / (1 + 2 * k_sol(3) * beta_d^2 / b_i^2);
            phi_i = 1 / (1 + 2 * k_sol(3) * (beta_d^2 + beta_c^2) / b_i^2);

            if Range.T_ritorno(t) == 50         % SLU spostamento
                im_c = Spost_limite(i,j);
            elseif Range.T_ritorno(t) == 1
                im_c = Acc_limite(i,j);
            end

            H_i_Cap = k_sol(1) * exp(-k_sol(3) * (log(im_c))^2 - k_sol(2) * log(im_c));

            FieldName3 = sprintf('Vel_%02d', Range.Vento(v));

            % if Range.T_ritorno(t) == 50         % SLU spostamento
            %     lambda_PS_spost.(FieldName3)(i,j) = sqrt(phi_i) * k_sol(1)^(1 - phi_i) * H_i_Cap^(phi_i) * ...
            %     exp(0.5 * q_i * (k_sol(2))^2 * (beta_c^2 + phi_i * beta_d^2));
            % elseif Range.T_ritorno(t) == 1
            %     lambda_PS_acc.(FieldName3)(i,j) = sqrt(phi_i) * k_sol(1)^(1 - phi_i) * H_i_Cap^(phi_i) * ...
            %     exp(0.5 * q_i * (k_sol(2))^2 * (beta_c^2 + phi_i * beta_d^2));
            % end

            MAF = sqrt(phi_i) * k_sol(1)^(1 - phi_i) * H_i_Cap^(phi_i) * ...
                exp(0.5 * q_i * (k_sol(2))^2 * (beta_c^2 + phi_i * beta_d^2));

        end
    end

end