classdef lossBasedDesignISO
%lossBasedDesignIsolatedStructures

properties
Selected
SelectedProp

EALtarget
ssType
isolationType
parameters
IMdef

Dyiso
fyiso
alphaiso
k1
k2
t1
t2
muiso 
DuIso
DncIso

hiso
aratio
sigma
areaiso
numiso

Riso
diamFPS
Teff
Deff

fySS
kSS
heffstar
expDisp
expVss
acc2DOF
accShape
PSDMSS
PFAGivenIM

powerLawDuct
powerLawAccRat

ductThresholdsIso

fragIsoMedian
fragIsoStDev
fragCollapse

fragSSMedian
fragSSStDev

vulnDisaggregated
labelsVuln
vulnerabilitiesIso
vulnerabilitiesSS
vulnerabilitiesNSCD
vulnerabilitiesNSCA
vulnerabilitiesCollapse

ealIso
ealSS
ealNSCD
ealNSCA
ealCollapse
eal

mafeDSColIso
mafeDSyieldSS

det = struct('Alead',[],'Arubber',[],'HeightLRB',[])

isEALTARGET
isISOCOLACC
isSSYIELDACC
isCANDIDATE
indCANDIDATES
designpar

end

properties(Access = 'private')

hyst
isExtrapolated    

end

methods % Macro-methods

function self = lossBasedDesignISO...
        (isolationtype, sstype,options)
    
    % Set some variables
    % self.EALtarget=EALtarget;
    % options.General.EALtarget=EALtarget;
    self.isolationType = isolationtype;
    self.ssType= sstype;
    
    % Set all parameters and define IM sampling
    self = setAllParameters(self, options);
    fragV=self.parameters.FragVuln;
    self.IMdef = linspace(0, fragV.maxIM,fragV.samplesIM)';
    
    % Set seed structures (design space)
    DET=self.parameters.Controls.DET;
    if   DET==0, self = setSeedSDoFproperties(self); 
        elseif   DET==1, self = setSeedSDoFproperties2(self); 
            else, warning('Specify 1 or 0 for parameters.Control.DET'); 
    end

    self.EALtarget=self.parameters.General.EALtarget;

end

function self = getSSparameters(self,manualDeltaY)
    
    if nargin < 2; manualDeltaY = 0; end
    
    SS = self.parameters.SS;
    SS.StoreyHeights=cumsum(SS.InterStoreyHeights);
    SS.FloorsMass=SS.FloorsWeight./9.81;
    SS.Gamma=sum(SS.FloorsWeight(2:end))/sum(SS.FloorsWeight);
    SS.numStoreys=length(SS.FloorsWeight)-1;
    self.parameters.SS=SS; 
    
    self = DBDprop(self,manualDeltaY);
    
end

function self = getSeedEAL(self, GPfullFit,GPfullFitAcc)
    


    self = setLossAssessmentMode(self);
    self = getPowerLawDuct(self, GPfullFit); %PL mod needed
    self = getPowerLawAcc(self, GPfullFitAcc);  %PL mod needed 
    self = calcAccProfile(self);
%     self = calcHeffstar(self);   
%     self = calcSSPSDM(self);
%     self = ductSSdonatello(self);
    self = ductSSYe(self); %PL mod needed
    self = setIsoDSThreshold(self);     
    
    self = getFragIso(self,GPfullFit); %PL mod needed
    self = getfragCollapse(self);
    self = getFragSS(self); 
        
    self = getVulnISO(self);
    self = getVulnSS(self);
    self = getVulnNSCD(self);
    self = getVulnNSCA(self,GPfullFitAcc.PSDMformula);  
    self = getVulnCollapse(self);
    
    self = getEALs(self); 

end

function self= getSeedMAFEds(self)

    haz = self.parameters.Hazard;

    for n = numel(self.fyiso): -1 : 1    
        for l = size(haz.intensityHazard,1) : -1 : 1
            intHaz(l) = interp1(haz.periodsHazard, ...
                haz.intensityHazard(l,:), self.t1(n), ...
                'linear', 'extrap');
        end

        hazCurve = self.logExtrapHazardCurve(...
            [intHaz', haz.MAFEhazard(:)], haz.faultRate, ...
            haz.extrapMode);

        self.mafeDSColIso(n,:) = self.calculateMAFEds(...
            self.fragIsoMedian(n,end), self.fragIsoStDev(n,end), hazCurve);

        for SS=length(self.fySS):-1:1
            self.mafeDSyieldSS(n,SS) = self.calculateMAFEds(...
                self.fragSSMedian(n,end,SS), self.fragSSStDev(n,end,SS), hazCurve);
        end
    end

    self.mafeDSColIso(isnan(self.mafeDSColIso))=0;
    self.mafeDSyieldSS(isnan(self.mafeDSyieldSS))=0;
end

function self = getDesignCandidates(self)

    % To be included a check for effectiveness of the isolation
    % system (check DuIso/DySS T1SS/T1Iso T1SS/TeffIso)
    Gen=self.parameters.General;
    self.isISOCOLACC = ...
         self.mafeDSColIso(:,end)...
         <= Gen.MAFECollapseLim;   

     for SS=length(self.fySS):-1:1
         self.isSSYIELDACC(:,SS) = self.mafeDSyieldSS(:,SS)...
         <= Gen.MAFESSYieldLim; 
         self.isEALTARGET(:,SS) = abs(self.eal(:,SS) - Gen.EALtarget) <= ...
            Gen.toleranceEAL;
         self.isCANDIDATE(:,SS) = self.isEALTARGET(:,SS) & ...
             self.isISOCOLACC & self.isSSYIELDACC(:,SS);
         self.indCANDIDATES{1,SS} = find(self.isCANDIDATE(:,SS)==1);
     end  
     
    self=getLRBDetailing(self);
    self=candidatesInfo(self);
    
end

function self = selectDesign(self,SelectedDsgn)
    self.Selected=SelectedDsgn;
    if SelectedDsgn==0
        self.Selected = input('Provide Select Design index [iso SS]: ');
    end
    switch self.isolationType
        case 'LRB'
%             printDesigninfoLRB(self)

        case 'FPS'
            printDesigninfoFPS(self)
    end
end

function self = getDesignProperties(self)
    
    sel=self.Selected;
    SP.t1SDoF=self.t1(sel(1));
    SP.t2=self.t2(sel(1));
    SP.dyIso=self.Dyiso(sel(1));
    SP.Fyiso=self.fyiso(sel(1));
    SP.FySS=self.fySS(sel(2));
    SP.k1iso=self.k1(sel(1));
    SP.Alphaiso=self.alphaiso(sel(1));
    SP.interSH=self.parameters.SS.InterStoreyHeights;
    SP.FloorsWeight=self.parameters.SS.FloorsWeight;
    SP.Heff=self.parameters.SS.heff;
    SP.Drifty=self.parameters.SS.drifty; 
    SP.dySStop=self.parameters.SS.DYatTop;
    SP.numStoreys=self.parameters.SS.numStoreys;
    self.SelectedProp=SP;


end

function self = getLRBDetailing(self)
        
    par=self.parameters;
    Wtot=sum(par.SS.FloorsWeight);
    self.det.Alead = (self.fyiso.*Wtot)./...
        (par.Iso.fylead+(par.Iso.Glead.*self.alphaiso.*par.Iso.epsYlead)./(1-self.alphaiso))./1000;
    self.det.Arubber = self.alphaiso.*self.det.Alead*par.Iso.Glead./(par.Iso.Grubber.*(1-self.alphaiso));
    self.det.HeightLRB = (par.Iso.Glead*1000.*self.det.Alead+par.Iso.Grubber*1000.*self.det.Arubber).*self.t1.^2*9.81./...
                         (4*pi^2*Wtot);
end

function plotCandidates(self)

    plotCandidatesScatter(self,ones(length(self.fyiso),length(self.fySS)),'All seeds ');
    plotCandidatesScatter(self,self.isISOCOLACC.*self.isSSYIELDACC.*ones(length(self.fyiso),length(self.fySS)),'Collapse OK ');
    plotCandidatesScatter(self,self.isCANDIDATE,'Design Candidates ');
%     plotCandidatesBackbone(self);

end 



end

methods % Micro-methods
    
    function self = setLossAssessmentMode(self)

            storeyLabels = {'(s1)', '(s2)', '(s3)', '(s4)', '(s5)', ...
                '(s6)', '(s7)', '(s8)', '(s9)', '(s10)', '(s11)', '(s12)'};
            
            self.labelsVuln = {'SS'};
            self.labelsVuln(2) = {'Iso'};

            if self.parameters.NSC.explicitNSC
                self = buildStoreyLossFx(self);
                                
                self.labelsVuln(end+1: ...
                    end+self.parameters.SS.numStoreys) = ...
                    strcat('NSCD', ...
                    storeyLabels(1:self.parameters.SS.numStoreys));
                self.labelsVuln(end+1: ...
                    end+self.parameters.SS.numStoreys) = ...
                    strcat('NSCA', ...
                    storeyLabels(1:self.parameters.SS.numStoreys));
            end

            self.labelsVuln(end+1) = {'Indirect'};

            self = getNormDLRs(self); 
    end
        
       function self = buildStoreyLossFx(self)

        num = [self.parameters.SS.numStoreys,...
               self.parameters.SS.numStoreys+1];
            groups = {'Dsens', 'Asens'};

            for g = 1 : 2
                for n = num(g) : -1 : 1
                    p = self.parameters.NSC.(groups{g})(n).SLFpars;
                    switch self.parameters.General.lossType
                        case 'directLosses'
                            factor=self.parameters.NSC.(groups{g})(n).costRatio.*100;
            
                        case 'repairTime'
                            F=self.parameters.LossParms;
                            factor=self.parameters.NSC.(groups{g})(n).totalRepairTime/...
                                  (F.numWorkersPerFloor*...
                                   F.workingHoursRatio*...
                                   F.workDaysRatio);
                    end
                    switch self.parameters.NSC.(groups{g})(n).SLFshape
                        case 'Papadopoulos'
                            self.parameters.NSC.(groups{g})(n).SLF = @(edp) factor.* ...
                                ( p(5) * (edp.^p(1) ./ (p(2)^p(1) + edp.^p(1))) + ...
                                (1-p(5)) * (edp.^p(3) ./ (p(4)^p(3) + edp.^p(3))) )...
                                .* (( p(5) * (edp.^p(1) ./ (p(2)^p(1) + edp.^p(1))) + ...
                                (1-p(5)) * (edp.^p(3) ./ (p(4)^p(3) + edp.^p(3))) )>0);
                        case 'Weibull'
                            self.parameters.NSC.(groups{g})(n).SLF = @(edp) factor.*...
                                (p(1).*(1-exp(-(edp/p(2)).^p(3))));
                    end
                end
            end
       end

     function self = getNormDLRs(self)

        switch self.parameters.General.lossType
            case 'directLosses'
                factorIso=self.parameters.Iso.costRatio*100;
                factorSS=self.parameters.SS.costRatio*100;
                factorCollapse=100;

            case 'repairTime'
                F=self.parameters.LossParms;
                factorIso=1/(F.numWorkersPerFloor*...
                       F.workingHoursRatio*...
                       F.workDaysRatio);
                factorSS=1/(F.numWorkersPerFloor*...
                       F.workingHoursRatio*...
                       F.workDaysRatio);
                factorCollapse=factorIso;
        end

        self.parameters.Iso.DLRnorm=self.parameters.Iso.DLR*factorIso;
        self.parameters.Iso.DLRcollapse=self.parameters.Iso.DLR(end)*factorCollapse;
        self.parameters.SS.DLRnorm=self.parameters.SS.DLR*factorSS;
    end
        
    function self = DBDprop(self,manualDeltaY)
        
        SS=self.parameters.SS;
        Htot=SS.StoreyHeights(end)-SS.StoreyHeights(1);
        
        switch self.ssType
            case 'Frame'

                if manualDeltaY == 0
                    if strcmp(SS.Material,'Steel'), k=0.65;
                    elseif strcmp(SS.Material,'Concrete'), k=0.5; end
                    SS.drifty=k*SS.epsYsteel*SS.lengthBeams/SS.heightBeams;
                else 
                    SS.drifty=manualDeltaY/SS.StoreyHeights(end);
                end 

                %Assuming a linear distribution of deformations  
                SS.dY_func=@(h) SS.drifty*h; 
                SS.dispShape=@(h) h/SS.storeyHeights;

                SS.heff=sum(SS.FloorsWeight(2:end).*SS.StoreyHeights(2:end).^2)/...
                    sum(SS.FloorsWeight(2:end).*SS.StoreyHeights(2:end))-SS.StoreyHeights(1);
%                 SS.effDispGivenUnitRoofDisp = ...
%                 sum(SS.FloorsMass .* (frame.dispShape.^2)) / ...
%                 sum(SS.FloorsMass .* frame.dispShape);
                                    
            case 'Wall'
                if manualDeltaY == 0
                    curvY=2*SS.epsYsteel./SS.lengthWall*0.5; % half of yield dip. reported by "DDBD of structures..." (Priestly...)
                else 
                    curvY=3/2*manualDeltaY/((Htot))^2; 
                end 
                
                %Assuming a cantilever displacement profile
                SS.dY_func=@(h) curvY.*h.^2.*(1-h./(3*(Htot))); 
                SS.heff=sum(SS.FloorsWeight(2:end).*...
                    SS.dY_func(SS.StoreyHeights(2:end)-SS.StoreyHeights(1)).*...
                    (SS.StoreyHeights(2:end)-SS.StoreyHeights(1)))./...
                    sum(SS.FloorsWeight(2:end).*...
                    SS.dY_func(SS.StoreyHeights(2:end)-SS.StoreyHeights(1)));
        end
        
        SS.Dy=SS.dY_func(SS.heff);
        SS.drifty=diff([0 SS.dY_func(SS.StoreyHeights(2:end)-...
            SS.StoreyHeights(1))./SS.InterStoreyHeights(2:end)]); 
        SS.kSS=self.fySS/SS.Dy;
        self.parameters.SS=SS;
    end

    function self = calcAccProfile(self)

    switch self.parameters.SS.accProf
        case 'FEMA'
            self = calcAccShapeFEMA(self);
        case 'Linear'
            self = calcAcc2DOFModal(self);
            self = calcAccShapeLinear(self);
        case 'LinearElasticProp'
            self=calcAcc2DOFModalLinearEalstic(self);
            self = calcAccShapeLinear(self);
        case 'Miranda'
            self = calcAccShapeMiranda(self);
        case 'Constant'
            self = calcAccShapeConst(self);
    end

    end
      
    function self = calcAccShapeFEMA(self)
        % FEMA P58 - Volume 1 - Eq. 5-11

        switch self.ssType
            case 'Frame'
                a0 = 0.66; a1 = -0.25; a2 = -0.080;
                a3 = -0.039; a4 = 0; a5 = 0;
            case 'Wall'
                a0 = 0.66; a1 = -0.15; a2 = -0.084;
                a3 = -0.26; a4 = 0.57; a5 = 0;
        end

        heights = self.parameters.SS.StoreyHeights';

        hRatio = (heights(2:end)-heights(1)) ./ ...
            (heights(end)-heights(1));

        self.parameters.(self.ssType).accShape = ...
            @(Tsdof,SAoverfy)exp(a0 + a1*Tsdof + a2*SAoverfy + ...
            a3*hRatio + a4*hRatio.^2 + a5*hRatio.^3);

        for SSt=1:size(self.fySS)
                AccShape(:,SSt) = ...
                    self.parameters.(self.ssType).accShape...
                    (self.parameters.SS.t1(SSt),1);
        end

        self.accShape=[1;AccShape].*ones(length(self.IMdef),1);

    end

      function self = calcAcc2DOFModal(self)
        % acc profile using the substitute structure and the effective properties
        % to calc acc at heff. Then computtion of floor acc using the SS
        % disp profile.

        SS=self.parameters.SS;
        M=[SS.FloorsMass(1),0 ; ...
            0,sum(SS.FloorsMass(2:end))*SS.meffFactor];

    for SSt=numel(self.fySS):-1:1
        for n = numel(self.fyiso):-1:1 
            R=self.IMdef./self.fyiso(n);
            EC = R <= 1; %ElasticCases
            Q=self.fyiso(n)-(self.k2(n)*self.Dyiso(n));
            Diso(EC)=R(EC).*self.Dyiso(n);
            Diso(~EC)=self.Dyiso(n).*(1+self.powerLawDuct(n,1)*(R(~EC)-1));
            Teff(EC)=self.t1(n);
            keff(~EC)=self.k2(n)+Q./Diso(~EC);
            Teff(~EC)=2*pi.*(1./(keff(~EC)*9.81)).^0.5;
            Trat=Teff./SS.t1(SSt);
            epsilonEff = 1./(Trat.^2).*sum(SS.FloorsWeight)/...
                (SS.FloorsWeight(1)+sum(SS.FloorsWeight(2:end))*SS.meffFactor);
            modalShape = [ones(length(self.IMdef),1),1+epsilonEff'];

            for m=1:length(self.IMdef)
                M1(m)=modalShape(m,:)*M*modalShape(m,:)';
                NmodalShape(m,:)=modalShape(m,:)./(M1(m)^0.5);
                MPF(m) = (NmodalShape(m,:)*M*ones(length(SS.numStoreys)+1,1));
            end

            self.acc2DOF(:,n,:,SSt)=NmodalShape.*MPF';
        end
        
    end
    end


    function self = calcAccShapeLinear(self)

        SS=self.parameters.SS;
    
    for SSt=numel(self.fySS):-1:1
        self.accShape(:,:,1,SSt) = self.acc2DOF(:,:,1,SSt);
        for st = 1 : self.parameters.SS.numStoreys 
            self.accShape(:,:,st+1,SSt) = self.acc2DOF(:,:,1,SSt) +...
                ((self.acc2DOF(:,:,2,SSt)-self.acc2DOF(:,:,1,SSt))./...
                SS.heff).* (SS.StoreyHeights(st+1)-SS.StoreyHeights(1));
        end
    end
    end

    function self = calcAcc2DOFModalLinearEalstic(self)

        SS=self.parameters.SS;
        M=[sum(SS.FloorsMass(2:end))*SS.meffFactor,0 ; ...
            0,SS.FloorsMass(1)];

        for SSt=numel(self.fySS):-1:1
            for n = numel(self.fyiso):-1:1
            Teff=self.fyiso(n)/self.Dyiso(n);
            Trat=Teff./SS.t1(SSt);
            epsilonEff = 1./(Trat.^2).*sum(SS.FloorsWeight)/...
                (SS.FloorsWeight(1)+sum(SS.FloorsWeight(2:end))*SS.meffFactor);
            modalShape = [1,1+epsilonEff];
            MPF = modalShape*M*ones(length(SS.numStoreys)+1,1)/...
                  (modalShape*M*modalShape');
            self.acc2DOF(:,n,:,SSt)=ones(length(self.IMdef),1).*modalShape*MPF;
            end
        end
    end

    function self = calcAccShapeConst(self)
        
        self.accShape=ones(length(self.IMdef),length(self.fyiso),self.parameters.SS.numStoreys+1,numel(self.fySS));

    end

    function self = getPowerLawDuct(self, GPfullFit)

        hystdummy = cell(numel(self.fyiso),1);
        hystdummy(:) = {self.isolationType};

        [self.powerLawDuct, self.isExtrapolated] = ...
            GPfullFit.getPowerLaw(hystdummy, self.t1, self.fyiso, self.alphaiso);

    end

    function self = getPowerLawAcc(self, GPfullFitAcc)

        hystdummy = cell(numel(self.fyiso),1);
        hystdummy(:) = {self.isolationType};

        [self.powerLawAccRat, self.isExtrapolated] = ...
            GPfullFitAcc.getPowerLawAcc(hystdummy, self.t1, self.fyiso, self.alphaiso);

    end

    function self = getFragIso(self, GPfullFit)

        [self.fragIsoMedian, dummyStDev] = ...
            GPfullFit.predictFragGPIso(self.fyiso,...
            self.ductThresholdsIso, self.powerLawDuct);

        self.fragIsoStDev=dummyStDev;
        % SDoF approximation for isolator dosent introduce much more dispersion             
        %             if isnan(self.parameters.FragVuln.fixedBeta)
        %                 self.fragIsoStDev = ( (dummyStDev(:,end) * ...
        %                     ones(1,size(self.ductThresholdsIso,2))).^2 + ...
        %                     self.parameters.FragVuln.betaSDoFtoMDoF.^2 ) .^ 0.5;
        %             else
        %                 self.fragIsoStDev = ones(size(dummyStDev)) * ...
        % %                     self.parameters.FragVuln.fixedBeta;
        %             end

        if any(self.isExtrapolated==1)
            warning('Extrapolation: some Seed SDoF exceeds the limits of the GP')
        end

    end


    function self = getfragCollapse(self)

        switch self.parameters.Controls.ColLosses

            case 1
              for n = numel(self.fyiso): -1 : 1
                self.fragCollapse(:,n) = logncdf(self.IMdef,...
                log(self.fragIsoMedian(n,end)), self.fragIsoStDev(n,end));
              end

            case 0 
              for n = numel(self.fyiso): -1 : 1
                self.fragCollapse(:,n) = zeros(length(self.IMdef),1);
              end
        end
    end

    function self = getFragSS(self) 

        % Get fragilities median and st.Dev
        for SS=length(self.fySS):-1:1
            for n = numel(self.fyiso): -1 : 1

                self.fragSSMedian(n,:,SS)=...
                    interp1(self.PSDMSS(2:end,n,SS),self.IMdef(2:end),...
                    self.parameters.SS.DS);

                % To include reference to frag formula for flexibility in changing it
                EC=self.fragSSMedian(n,:,SS)/self.fyiso(n)<=1;
                fragStDev(n,EC,SS) = 0.01;
                fragStDev(n,~EC,SS) = self.powerLawAccRat(n,2);    
            end

            % Get total disperssion
            if isnan(self.parameters.FragVuln.fixedBeta)
            self.fragSSStDev(:,:,SS) = ((fragStDev(:,:,SS)).^2 + ...
                (self.parameters.FragVuln.betaSDoFtoMDoF.*...
                ones(numel(self.fyiso),size(self.parameters.SS.DS,2))).^2 ) .^ 0.5;
            else
            self.fragSSStDev(:,:,SS) = ones(self.powerLawAccRat(:,end),...
                length(self.fySS)) * self.parameters.FragVuln.fixedBeta;
            end   
        end

        % Control plots Fragilities
        % figure()
        % hold on
        % SS=1;
        % n=544;
        % plot(self.IMdef,logncdf(self.IMdef,...
        %     log(self.fragSSMedian(n,1)), self.fragSSStDev(n,1,SS)));
        % plot(self.IMdef,logncdf(self.IMdef,...
        %     log(self.fragSSMedian(n,2)), self.fragSSStDev(n,2,SS)));
        % plot(self.IMdef,logncdf(self.IMdef,...
        %     log(self.fragSSMedian(n,3)), self.fragSSStDev(n,3,SS)));
        % title(strcat('SS Yield Fragility for fySS=',num2str(self.fySS(SS)),'g'));
        % xlabel('IM: Sa(t1) (g)')
        % ylabel('PE')

        % Control plots Vulneravilities
        % figure()
        % hold on
        % SS=1;
        % n=544;
        % plot(self.IMdef,...
        %     self.vulnerabilitiesSS(:,n,SS));
        % title(strcat('SS vulnerability fySS=',num2str(self.fySS(SS)),'g'));
        % xlabel('IM: Sa(t1) (g)')
        % ylabel('Loss Ratio')
    end
 
    function self = getVulnISO(self)
       
        for n = numel(self.fyiso): -1 : 1
            for ds = size(self.fragIsoMedian,2)-1 : -1 : 1
                fragilities(:,ds) = logncdf(self.IMdef,...
                    log(self.fragIsoMedian(n,ds)), self.fragIsoStDev(n,ds));
            end

        vulnerability = ...
            VULNERABILITYbuilding([self.IMdef, fragilities], ...
            self.parameters.Iso.DLRnorm(1:end-1), 'NOplot', ...
            'S_a(T_1)', 0, self.parameters.FragVuln.fxIndLosses, ...
            self.parameters.FragVuln.parIndLosses);
            self.vulnerabilitiesIso(:,n) = ...
            vulnerability(:,2) + vulnerability(:,3);
        end
        
    end
    
    function self = getVulnSS(self)
    
    % Compute the vulnerability curves
        for SS=length(self.fySS):-1:1
            for n = numel(self.fyiso): -1 : 1
                for ds = size(self.fragSSMedian,2) : -1 : 1
                    fragilities(:,ds) = logncdf(self.IMdef,...
                        log(self.fragSSMedian(n,ds,SS)), self.fragSSStDev(n,ds,SS));
                end
                vulnerability = ...
                    VULNERABILITYbuilding([self.IMdef, fragilities(:,:)], ...
                    self.parameters.SS.DLRnorm, 'noplot', ...
                    'S_a(T_1)', 0, self.parameters.FragVuln.fxIndLosses, ...
                    self.parameters.FragVuln.parIndLosses);
                    self.vulnerabilitiesSS(:,n,SS) = ...
                    vulnerability(:,2) + vulnerability(:,3);
            end
        end  
    end
    
    function self = getVulnNSCD(self)
            
         for SSt = 1 : numel(self.fySS)
            for n = 1 : numel(self.fyiso)
                driftGivenIM = ...
                    self.PSDMSS(:,n,SSt).*self.parameters.SS.drifty;
                for st = 1 : self.parameters.SS.numStoreys                    
                    self.vulnerabilitiesNSCD(:,n,SSt,st) = ...
                        self.parameters.NSC.Dsens(st).SLF(driftGivenIM(:,st).*100);
                    self.vulnerabilitiesNSCD(:,n,SSt,st) = ...
                        self.vulnerabilitiesNSCD(:,n,SSt,st).*(self.vulnerabilitiesNSCD(:,n,SSt,st)>0);
                end
            end
         end
    end
    
         
    function self = getVulnNSCA(self,PSDMformula)
         self.PFAGivenIM=[];   

         for SSt = 1 : numel(self.fySS)
            for n = 1 : numel(self.fyiso)
                R = self.IMdef./self.fyiso(n);
                EC = R <= 0.95;
                for st = self.parameters.SS.numStoreys+1 : -1 :1
                    self.PFAGivenIM(EC,n,SSt,st) = R(EC).*self.accShape(EC,n,st,SSt)*self.fyiso(n);
                    self.PFAGivenIM(~EC,n,SSt,st) = PSDMformula.AccRatio...
                        (R(~EC), self.powerLawAccRat(n,1), 0).*...
                        self.accShape(~EC,n,st,SSt).*self.fyiso(n); % mean
                    self.vulnerabilitiesNSCA(:,n,SSt,st) = ...
                        self.parameters.NSC.Asens(st).SLF(self.PFAGivenIM(:,n,SSt,st));
                    self.vulnerabilitiesNSCA(:,n,SSt,st) = ...
                        self.vulnerabilitiesNSCA(:,n,SSt,st).*(self.vulnerabilitiesNSCA(:,n,SSt,st)>0);
                end
            end
         end
    end

    function self = getVulnCollapse(self)

        for n = 1 : numel(self.fyiso)
            self.vulnerabilitiesCollapse(:,n) = ...
                    self.fragCollapse(:,n).*self.parameters.Iso.DLRcollapse;
        end
    end
    
    
    function self = getEALs(self)

        haz = self.parameters.Hazard;

        for n = numel(self.t1): -1 : 1

            for l = size(haz.intensityHazard,1) : -1 : 1
                intHaz(l) = interp1(haz.periodsHazard, ...
                    haz.intensityHazard(l,:), unique(self.t1(n)), ...
                    'linear', 'extrap');
            end

            hazCurve = self.logExtrapHazardCurve(...
                [intHaz', haz.MAFEhazard(:)], haz.faultRate, ...
                haz.extrapMode);

            self.ealIso(n,1) = EALcalculator(...
                [self.IMdef self.vulnerabilitiesIso(:,n)],hazCurve,...
                'Noplot',self.fragCollapse(:,n),0); 

            for SS=length(self.fySS):-1:1
                self.ealSS(n,SS) = EALcalculator(...
                    [self.IMdef self.vulnerabilitiesSS(:,n,SS)],hazCurve,...
                    'Noplot',self.fragCollapse(:,n),0); 

                for St= self.parameters.SS.numStoreys:-1:1
                    eal_NSCD(n,St,SS) = EALcalculator(...
                    [self.IMdef(2:end) self.vulnerabilitiesNSCD(2:end,n,SS,St)],...
                    hazCurve, 'NOplot',self.fragCollapse(:,n),0); 
                end

                for St= self.parameters.SS.numStoreys+1:-1:1 
                    eal_NSCA(n,St,SS) = EALcalculator(...
                    [self.IMdef self.vulnerabilitiesNSCA(:,n,SS,St)],...
                    hazCurve, 'NOplot',self.fragCollapse(:,n),0); 
                end                    
            end

            self.ealCollapse(n,1) = EALcalculator(...
            [self.IMdef self.vulnerabilitiesCollapse(:,n)],hazCurve,...
            'Noplot',zeros(length(self.IMdef),1),0);

        end

        self.ealNSCD = squeeze(sum(eal_NSCD,2));
        self.ealNSCA = squeeze(sum(eal_NSCA,2));

        self.ealIso(isnan(self.ealIso))=0;
        self.ealSS(isnan(self.ealSS))=0;
        self.ealNSCD(isnan(self.ealNSCD))=0;
        self.ealNSCA(isnan(self.ealNSCA))=0;

        self.eal=self.ealIso+self.ealSS+...
                 self.ealNSCD+self.ealNSCA+...
                 self.ealCollapse; 
    end
    
end

methods(Access = 'private')

    function self = setAllParameters(self, options)
    % setAllParameters deals with the optional parameters

    % build basic parameters
    macroFieldsPar = {'General', 'FragVuln', 'Hazard', 'Iso', 'SS','NSC','Controls','Loss','LossParms'};

    microFieldsPar{1} = {'EALtarget','toleranceEAL','MAFECollapseLim','MAFESSYieldLim', 'hysteresis','lossType'};
    microFieldsParVals{1} = {0.01, 0.005,4.1E-5,1E-4,'LRB','directLosses'};

    microFieldsPar{2} = {'fixedBeta','betaSDoFtoMDoF', 'maxIM', 'samplesIM', 'fxIndLosses', 'parIndLosses'};
    microFieldsParVals{2} = { NaN ,0, 2.5, 1000, @(Ldir,par)par(1)*normcdf(-0.2+0.4*Ldir,0,par(2)), [0, 0.15]};

    microFieldsPar{3} = {'faultRate', 'periodsHazard', 'extrapMode', 'MAFEhazard', 'intensityHazard', 'codeSpectra' };
    microFieldsParVals{3} = { 0.58, [0,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1,1.5,2], 'powerLaw', ...
        [0.033214624,0.019885045,0.013862944,0.009885926,0.007133499,0.004969227,0.00210721,0.001025866,0.000404054], ...
        [0.118500000000000,0.231901076160571,0.284400000000000,0.284400000000000,0.284400000000000,0.284400000000000,0.249647374490702,0.166427989012818,0.124820042771853,0.0832129082285958,0.0597891517804315;0.156000000000000,0.294972075787921,0.363792000000000,0.363792000000000,0.363792000000000,0.363792000000000,0.326382422300220,0.217582434028613,0.163185282252471,0.108789450714589,0.0815918939267367;0.183000000000000,0.341305827038763,0.420458740558145,0.424194000000000,0.424194000000000,0.424194000000000,0.387811782857298,0.258529430931683,0.193893977166246,0.129261174159895,0.0969454922187405;0.213000000000000,0.392400522703245,0.482100784054867,0.490752000000000,0.490752000000000,0.490752000000000,0.455921264959425,0.303930770035728,0.227943687239298,0.151960368720452,0.113969728298491;0.241667382400000,0.438980114035280,0.537636479852921,0.556076646902400,0.556076646902400,0.556076646902400,0.531683535429199,0.354443451456181,0.265829361119539,0.177218032623064,0.132913118895613;0.274027891000000,0.495860969580310,0.606777508870465,0.634374567665000,0.634374567665000,0.634374567665000,0.618336644289972,0.412206673094016,0.309150327439777,0.206097985927785,0.154572902412411;0.347077173600000,0.621962545909488,0.759405232064231,0.820490438390400,0.820490438390400,0.820490438390400,0.820490438390400,0.565248464301760,0.423929063051237,0.282615900631437,0.211961012200778;0.407159360000000,0.727702500307648,0.887974070461473,0.977182464000000,0.977182464000000,0.977182464000000,0.977182464000000,0.695119133694674,0.521330086813096,0.347548986783303,0.260660585294040;0.467092460800000,0.836572020006089,1.02131179960913,1.14811326864640,1.14811326864640,1.14811326864640,1.14811326864640,0.846507595616605,0.634870338466778,0.423241969970881,0.317430187426319], ...
        {[0,0.118500000000000,0;0.146294907964624,0.284400000000000,0.00151250661533860;0.438884723893872,0.284400000000000,0.0136125595380474;0.448884723893872,0.278064297650119,0.0139227220658583;0.458884723893872,0.272004730112316,0.0142328845936691;0.468884723893872,0.266203629836465,0.0145430471214799;0.478884723893872,0.260644804997119,0.0148532096492907;0.488884723893872,0.255313388565835,0.0151633721771015;0.498884723893872,0.250195705535313,0.0154735347049123;0.508884723893872,0.245279155798452,0.0157836972327231;0.518884723893872,0.240552110570414,0.0160938597605340;0.528884723893872,0.236003820561215,0.0164040222883448;0.538884723893872,0.231624334372483,0.0167141848161556;0.548884723893872,0.227404425814465,0.0170243473439664;0.558884723893872,0.223335529026053,0.0173345098717772;0.568884723893872,0.219409680437654,0.0176446723995880;0.578884723893872,0.215619466749481,0.0179548349273988;0.588884723893872,0.211957978210200,0.0182649974552097;0.598884723893872,0.208418766576414,0.0185751599830205;0.608884723893872,0.204995807214853,0.0188853225108313;0.618884723893872,0.201683464878706,0.0191954850386421;0.628884723893872,0.198476462749127,0.0195056475664529;0.638884723893872,0.195369854384171,0.0198158100942637;0.648884723893872,0.192358998261503,0.0201259726220745;0.658884723893872,0.189439534639328,0.0204361351498853;0.668884723893872,0.186607364492931,0.0207462976776962;0.678884723893872,0.183858630312810,0.0210564602055070;0.688884723893872,0.181189698575239,0.0213666227333178;0.698884723893872,0.178597143717755,0.0216767852611286;0.708884723893872,0.176077733470956,0.0219869477889394;0.718884723893872,0.173628415414547,0.0222971103167502;0.728884723893872,0.171246304640062,0.0226072728445610;0.738884723893872,0.168928672415408,0.0229174353723719;0.748884723893872,0.166672935757608,0.0232275979001827;0.758884723893872,0.164476647829945,0.0235377604279935;0.768884723893872,0.162337489088476,0.0238479229558043;0.778884723893872,0.160253259110554,0.0241580854836151;0.788884723893872,0.158221869044848,0.0244682480114259;0.798884723893872,0.156241334628397,0.0247784105392367;0.808884723893872,0.154309769721642,0.0250885730670476;0.818884723893872,0.152425380317137,0.0253987355948584;0.828884723893872,0.150586458981959,0.0257088981226692;0.838884723893872,0.148791379697609,0.0260190606504800;0.848884723893872,0.147038593064636,0.0263292231782908;0.858884723893872,0.145326621842258,0.0266393857061016;0.868884723893872,0.143654056795989,0.0269495482339124;0.878884723893872,0.142019552828739,0.0272597107617232;0.888884723893872,0.140421825373073,0.0275698732895341;0.898884723893872,0.138859647024276,0.0278800358173449;0.908884723893872,0.137331844395695,0.0281901983451557;0.918884723893872,0.135837295179404,0.0285003608729665;0.928884723893872,0.134374925396747,0.0288105234007773;0.938884723893872,0.132943706824573,0.0291206859285881;0.948884723893872,0.131542654584223,0.0294308484563990;0.958884723893872,0.130170824881377,0.0297410109842098;0.968884723893872,0.128827312885872,0.0300511735120206;0.978884723893872,0.127511250741461,0.0303613360398314;0.988884723893872,0.126221805696346,0.0306714985676422;0.998884723893872,0.124958178345992,0.0309816610954530;1.00888472389387,0.123719600980446,0.0312918236232638;1.01888472389387,0.122505336028984,0.0316019861510746;1.02888472389387,0.121314674595453,0.0319121486788855;1.03888472389387,0.120146935078206,0.0322223112066963;1.04888472389387,0.119001461868985,0.0325324737345071;1.05888472389387,0.117877624125520,0.0328426362623179;1.06888472389387,0.116774814613040,0.0331527987901287;1.07888472389387,0.115692448610196,0.0334629613179395;1.08888472389387,0.114629962875283,0.0337731238457503;1.09888472389387,0.113586814668899,0.0340832863735612;1.10888472389387,0.112562480829489,0.0343934489013720;1.11888472389387,0.111556456898465,0.0347036114291828;1.12888472389387,0.110568256291819,0.0350137739569936;1.13888472389387,0.109597409515389,0.0353239364848044;1.14888472389387,0.108643463421094,0.0356340990126152;1.15888472389387,0.107705980501688,0.0359442615404260;1.16888472389387,0.106784538221709,0.0362544240682368;1.17888472389387,0.105878728382483,0.0365645865960477;1.18888472389387,0.104988156519168,0.0368747491238585;1.19888472389387,0.104112441327984,0.0371849116516693;1.20888472389387,0.103251214121865,0.0374950741794801;1.21888472389387,0.102404118312902,0.0378052367072909;1.22888472389387,0.101570808920070,0.0381153992351017;1.23888472389387,0.100750952100778,0.0384255617629126;1.24888472389387,0.0999442247049409,0.0387357242907234;1.25888472389387,0.0991503138502933,0.0390458868185342;1.26888472389387,0.0983689165177915,0.0393560493463450;1.27888472389387,0.0975997391659948,0.0396662118741558;1.28888472389387,0.0968424973633988,0.0399763744019666;1.29888472389387,0.0960969154377520,0.0402865369297774;1.30888472389387,0.0953627261414489,0.0405966994575882;1.31888472389387,0.0946396703321443,0.0409068619853991;1.32888472389387,0.0939274966677889,0.0412170245132099;1.33888472389387,0.0932259613153306,0.0415271870410207;1.34888472389387,0.0925348276723740,0.0418373495688315;1.35888472389387,0.0918538661011289,0.0421475120966423;1.36888472389387,0.0911828536740208,0.0424576746244531;1.37888472389387,0.0905215739303702,0.0427678371522639;1.38888472389387,0.0898698166435841,0.0430779996800747;1.39888472389387,0.0892273775983322,0.0433881622078856;1.40888472389387,0.0885940583772129,0.0436983247356964;1.41888472389387,0.0879696661564405,0.0440084872635072;1.42888472389387,0.0873540135101115,0.0443186497913180;1.43888472389387,0.0867469182226328,0.0446288123191288;1.44888472389387,0.0861482031089176,0.0449389748469396;1.45888472389387,0.0855576958419761,0.0452491373747504;1.46888472389387,0.0849752287875488,0.0455592999025612;1.47888472389387,0.0844006388454483,0.0458694624303721;1.48888472389387,0.0838337672972957,0.0461796249581829;1.49888472389387,0.0832744596603514,0.0464897874859937;1.50888472389387,0.0827225655471587,0.0467999500138045;1.51888472389387,0.0821779385307312,0.0471101125416153;1.52888472389387,0.0816404360150318,0.0474202750694261;1.53888472389387,0.0811099191105008,0.0477304375972370;1.54888472389387,0.0805862525144058,0.0480406001250478;1.55888472389387,0.0800693043957975,0.0483507626528586;1.56888472389387,0.0795589462848646,0.0486609251806694;1.57888472389387,0.0790550529664933,0.0489710877084802;1.58888472389387,0.0785575023778468,0.0492812502362910;1.59888472389387,0.0780661755097876,0.0495914127641018;1.60888472389387,0.0775809563119767,0.0499015752919127;1.61888472389387,0.0771017316014898,0.0502117378197235;1.62888472389387,0.0766283909747991,0.0505219003475343;1.63888472389387,0.0761608267229783,0.0508320628753451;1.64888472389387,0.0756989337499927,0.0511422254031559;1.65888472389387,0.0752426094939449,0.0514523879309667;1.66888472389387,0.0747917538511514,0.0517625504587775;1.67888472389387,0.0743462691029330,0.0520727129865883;1.68888472389387,0.0739060598450062,0.0523828755143991;1.69888472389387,0.0734710329193675,0.0526930380422100;1.70888472389387,0.0730410973485704,0.0530032005700208;1.71888472389387,0.0726161642722958,0.0533133630978316;1.72888472389387,0.0721961468861236,0.0536235256256424;1.73888472389387,0.0717809603824176,0.0539336881534532;1.74888472389387,0.0713705218932381,0.0542438506812640;1.75888472389387,0.0709647504352017,0.0545540132090749;1.76888472389387,0.0705635668562119,0.0548641757368857;1.77888472389387,0.0701668937839864,0.0551743382646965;1.78888472389387,0.0697746555763100,0.0554845007925073;1.79888472389387,0.0693867782729479,0.0557946633203181;1.80888472389387,0.0690031895491536,0.0561048258481289;1.81888472389387,0.0686238186707099,0.0564149883759397;1.82888472389387,0.0682485964504455,0.0567251509037505;1.83888472389387,0.0678774552061704,0.0570353134315614;1.84888472389387,0.0675103287199759,0.0573454759593722;1.85888472389387,0.0671471521988490,0.0576556384871830;1.86888472389387,0.0667878622365503,0.0579658010149938;1.87888472389387,0.0664323967767102,0.0582759635428046;1.88888472389387,0.0660806950770969,0.0585861260706154;1.89888472389387,0.0657326976750134,0.0588962885984262;1.90888472389387,0.0653883463537826,0.0592064511262371;1.91700000000000,0.0650775711731428,0.0594271403285519;1.92700000000000,0.0644038948476958,0.0594271403285519;1.93700000000000,0.0637406253758758,0.0594271403285519;1.94700000000000,0.0630875495048417,0.0594271403285519;1.95700000000000,0.0624444594162588,0.0594271403285519;1.96700000000000,0.0618111525609515,0.0594271403285519;1.97700000000000,0.0611874314993952,0.0594271403285519;1.98700000000000,0.0605731037478131,0.0594271403285519;1.99700000000000,0.0599679816296551,0.0594271403285519;2.00700000000000,0.0593718821322430,0.0594271403285519;2.01700000000000,0.0587846267683784,0.0594271403285519;2.02700000000000,0.0582060414427185,0.0594271403285519;2.03700000000000,0.0576359563227323,0.0594271403285519;2.04700000000000,0.0570742057140585,0.0594271403285519;2.05700000000000,0.0565206279400951,0.0594271403285519;2.06700000000000,0.0559750652256564,0.0594271403285519;2.07700000000000,0.0554373635845419,0.0594271403285519;2.08700000000000,0.0549073727108673,0.0594271403285519;2.09700000000000,0.0543849458740134,0.0594271403285519;2.10700000000000,0.0538699398170582,0.0594271403285519;2.11700000000000,0.0533622146585583,0.0594271403285519;2.12700000000000,0.0528616337975551,0.0594271403285519;2.13700000000000,0.0523680638216865,0.0594271403285519;2.14700000000000,0.0518813744182857,0.0594271403285519;2.15700000000000,0.0514014382883599,0.0594271403285519;2.16700000000000,0.0509281310633406,0.0594271403285519;2.17700000000000,0.0504613312245044,0.0594271403285519;2.18700000000000,0.0500009200249676,0.0594271403285519;2.19700000000000,0.0495467814141598,0.0594271403285519;2.20700000000000,0.0490988019646882,0.0594271403285519;2.21700000000000,0.0486568708015052,0.0594271403285519;2.22700000000000,0.0482208795332983,0.0594271403285519;2.23700000000000,0.0477907221860212,0.0594271403285519;2.24700000000000,0.0473662951384914,0.0594271403285519;2.25700000000000,0.0469474970599810,0.0594271403285519;2.26700000000000,0.0465342288497299,0.0594271403285519;2.27700000000000,0.0461263935783142,0.0594271403285519;2.28700000000000,0.0457238964308062,0.0594271403285519;2.29700000000000,0.0453266446516617,0.0594271403285519;2.30700000000000,0.0449345474912766,0.0594271403285519;2.31700000000000,0.0445475161541543,0.0594271403285519;2.32700000000000,0.0441654637486289,0.0594271403285519;2.33700000000000,0.0437883052380917,0.0594271403285519;2.34700000000000,0.0434159573936684,0.0594271403285519;2.35700000000000,0.0430483387482991,0.0594271403285519;2.36700000000000,0.0426853695521739,0.0594271403285519;2.37700000000000,0.0423269717294772,0.0594271403285519;2.38700000000000,0.0419730688363989,0.0594271403285519;2.39700000000000,0.0416235860203678,0.0594271403285519;2.40700000000000,0.0412784499804699,0.0594271403285519;2.41700000000000,0.0409375889290090,0.0594271403285519;2.42700000000000,0.0406009325541747,0.0594271403285519;2.43700000000000,0.0402684119837802,0.0594271403285519;2.44700000000000,0.0399399597500354,0.0594271403285519;2.45700000000000,0.0396155097553209,0.0594271403285519;2.46700000000000,0.0392949972389328,0.0594271403285519;2.47700000000000,0.0389783587447634,0.0594271403285519;2.48700000000000,0.0386655320898911,0.0594271403285519;2.49700000000000,0.0383564563340485,0.0594271403285519;2.50700000000000,0.0380510717499417,0.0594271403285519;2.51700000000000,0.0377493197943929,0.0594271403285519;2.52700000000000,0.0374511430802810,0.0594271403285519;2.53700000000000,0.0371564853492551,0.0594271403285519;2.54700000000000,0.0368652914451961,0.0594271403285519;2.55700000000000,0.0365775072884039,0.0594271403285519;2.56700000000000,0.0362930798504860,0.0594271403285519;2.57700000000000,0.0360119571299286,0.0594271403285519;2.58700000000000,0.0357340881283255,0.0594271403285519;2.59700000000000,0.0354594228272484,0.0594271403285519;2.60700000000000,0.0351879121657353,0.0594271403285519;2.61700000000000,0.0349195080183812,0.0594271403285519;2.62700000000000,0.0346541631740110,0.0594271403285519;2.63700000000000,0.0343918313149171,0.0594271403285519;2.64700000000000,0.0341324669966455,0.0594271403285519;2.65700000000000,0.0338760256283137,0.0594271403285519;2.66700000000000,0.0336224634534434,0.0594271403285519;2.67700000000000,0.0333717375312938,0.0594271403285519;2.68700000000000,0.0331238057186810,0.0594271403285519;2.69700000000000,0.0328786266522670,0.0594271403285519;2.70700000000000,0.0326361597313072,0.0594271403285519;2.71700000000000,0.0323963651008406,0.0594271403285519;2.72700000000000,0.0321592036353115,0.0594271403285519;2.73700000000000,0.0319246369226084,0.0594271403285519;2.74700000000000,0.0316926272485097,0.0594271403285519;2.75700000000000,0.0314631375815232,0.0594271403285519;2.76700000000000,0.0312361315581086,0.0594271403285519;2.77700000000000,0.0310115734682714,0.0594271403285519;2.78700000000000,0.0307894282415190,0.0594271403285519;2.79700000000000,0.0305696614331663,0.0594271403285519;2.80700000000000,0.0303522392109831,0.0594271403285519;2.81700000000000,0.0301371283421727,0.0594271403285519;2.82700000000000,0.0299242961806717,0.0594271403285519;2.83700000000000,0.0297137106547635,0.0594271403285519;2.84700000000000,0.0295053402549951,0.0594271403285519;2.85700000000000,0.0292991540223895,0.0594271403285519;2.86700000000000,0.0290951215369462,0.0594271403285519;2.87700000000000,0.0288932129064195,0.0594271403285519;2.88700000000000,0.0286933987553704,0.0594271403285519;2.89700000000000,0.0284956502144803,0.0594271403285519;2.90700000000000,0.0282999389101239,0.0594271403285519;2.91700000000000,0.0281062369541898,0.0594271403285519;2.92700000000000,0.0279145169341459,0.0594271403285519;2.93700000000000,0.0277247519033397,0.0594271403285519;2.94700000000000,0.0275369153715297,0.0594271403285519;2.95700000000000,0.0273509812956399,0.0594271403285519;2.96700000000000,0.0271669240707324,0.0594271403285519;2.97700000000000,0.0269847185211918,0.0594271403285519;2.98700000000000,0.0268043398921159,0.0594271403285519;2.99700000000000,0.0266257638409068,0.0594271403285519;3.00700000000000,0.0264489664290582,0.0594271403285519;3.01700000000000,0.0262739241141321,0.0594271403285519;3.02700000000000,0.0261006137419211,0.0594271403285519;3.03700000000000,0.0259290125387913,0.0594271403285519;3.04700000000000,0.0257590981042003,0.0594271403285519;3.05700000000000,0.0255908484033865,0.0594271403285519;3.06700000000000,0.0254242417602252,0.0594271403285519;3.07700000000000,0.0252592568502467,0.0594271403285519;3.08700000000000,0.0250958726938122,0.0594271403285519;3.09700000000000,0.0249340686494445,0.0594271403285519;3.10700000000000,0.0247738244073076,0.0594271403285519;3.11700000000000,0.0246151199828339,0.0594271403285519;3.12700000000000,0.0244579357104922,0.0594271403285519;3.13700000000000,0.0243022522376960,0.0594271403285519;3.14700000000000,0.0241480505188462,0.0594271403285519;3.15700000000000,0.0239953118095058,0.0594271403285519;3.16700000000000,0.0238440176607039,0.0594271403285519;3.17700000000000,0.0236941499133635,0.0594271403285519;3.18700000000000,0.0235456906928533,0.0594271403285519;3.19700000000000,0.0233986224036570,0.0594271403285519;3.20700000000000,0.0232529277241600,0.0594271403285519;3.21700000000000,0.0231085896015484,0.0594271403285519;3.22700000000000,0.0229655912468193,0.0594271403285519;3.23700000000000,0.0228239161298982,0.0594271403285519;3.24700000000000,0.0226835479748618,0.0594271403285519;3.25700000000000,0.0225444707552632,0.0594271403285519;3.26700000000000,0.0224066686895576,0.0594271403285519;3.27700000000000,0.0222701262366244,0.0594271403285519;3.28700000000000,0.0221348280913859,0.0594271403285519;3.29700000000000,0.0220007591805180,0.0594271403285519;3.30700000000000,0.0218679046582516,0.0594271403285519;3.31700000000000,0.0217362499022630,0.0594271403285519;3.32700000000000,0.0216057805096500,0.0594271403285519;3.33700000000000,0.0214764822929928,0.0594271403285519;3.34700000000000,0.0213483412764968,0.0594271403285519;3.35700000000000,0.0212213436922160,0.0594271403285519;3.36700000000000,0.0210954759763542,0.0594271403285519;3.37700000000000,0.0209707247656440,0.0594271403285519;3.38700000000000,0.0208470768937990,0.0594271403285519;3.39700000000000,0.0207245193880399,0.0594271403285519;3.40700000000000,0.0206030394656919,0.0594271403285519;3.41700000000000,0.0204826245308515,0.0594271403285519;3.42700000000000,0.0203632621711210,0.0594271403285519;3.43700000000000,0.0202449401544099,0.0594271403285519;3.44700000000000,0.0201276464258009,0.0594271403285519;3.45700000000000,0.0200113691044795,0.0594271403285519;3.46700000000000,0.0198960964807248,0.0594271403285519;3.47700000000000,0.0197818170129622,0.0594271403285519;3.48700000000000,0.0196685193248732,0.0594271403285519;3.49700000000000,0.0195561922025652,0.0594271403285519;3.50700000000000,0.0194448245917956,0.0594271403285519;3.51700000000000,0.0193344055952528,0.0594271403285519;3.52700000000000,0.0192249244698899,0.0594271403285519;3.53700000000000,0.0191163706243117,0.0594271403285519;3.54700000000000,0.0190087336162128,0.0594271403285519;3.55700000000000,0.0189020031498668,0.0594271403285519;3.56700000000000,0.0187961690736636,0.0594271403285519;3.57700000000000,0.0186912213776958,0.0594271403285519;3.58700000000000,0.0185871501913913,0.0594271403285519;3.59700000000000,0.0184839457811930,0.0594271403285519;3.60700000000000,0.0183815985482822,0.0594271403285519;3.61700000000000,0.0182800990263469,0.0594271403285519;3.62700000000000,0.0181794378793929,0.0594271403285519;3.63700000000000,0.0180796058995965,0.0594271403285519;3.64700000000000,0.0179805940051993,0.0594271403285519;3.65700000000000,0.0178823932384422,0.0594271403285519;3.66700000000000,0.0177849947635397,0.0594271403285519;3.67700000000000,0.0176883898646919,0.0594271403285519;3.68700000000000,0.0175925699441347,0.0594271403285519;3.69700000000000,0.0174975265202272,0.0594271403285519;3.70700000000000,0.0174032512255738,0.0594271403285519;3.71700000000000,0.0173097358051833,0.0594271403285519;3.72700000000000,0.0172169721146617,0.0594271403285519;3.73700000000000,0.0171249521184383,0.0594271403285519;3.74700000000000,0.0170336678880262,0.0594271403285519;3.75700000000000,0.0169431116003139,0.0594271403285519;3.76700000000000,0.0168532755358893,0.0594271403285519;3.77700000000000,0.0167641520773947,0.0594271403285519;3.78700000000000,0.0166757337079117,0.0594271403285519;3.79700000000000,0.0165880130093765,0.0594271403285519;3.80700000000000,0.0165009826610237,0.0594271403285519;3.81700000000000,0.0164146354378592,0.0594271403285519;3.82700000000000,0.0163289642091601,0.0594271403285519;3.83700000000000,0.0162439619370029,0.0594271403285519;3.84700000000000,0.0161596216748182,0.0594271403285519;3.85700000000000,0.0160759365659708,0.0594271403285519;3.86700000000000,0.0159928998423666,0.0594271403285519;3.87700000000000,0.0159105048230841,0.0594271403285519;3.88700000000000,0.0158287449130303,0.0594271403285519;3.89700000000000,0.0157476136016210,0.0594271403285519;3.90700000000000,0.0156671044614848,0.0594271403285519;3.91700000000000,0.0155872111471900,0.0594271403285519;3.92700000000000,0.0155079273939943,0.0594271403285519;3.93700000000000,0.0154292470166166,0.0594271403285519;3.94700000000000,0.0153511639080304,0.0594271403285519;3.95700000000000,0.0152736720382793,0.0594271403285519;3.96700000000000,0.0151967654533122,0.0594271403285519;3.97700000000000,0.0151204382738399,0.0594271403285519;3.98700000000000,0.0150446846942115,0.0594271403285519;3.99700000000000,0.0149694989813100,0.0594271403285519;4.00700000000000,0.0148948754734679,0.0594271403285519;4.01700000000000,0.0148208085794013,0.0594271403285519;4.02700000000000,0.0147472927771624,0.0594271403285519;4.03700000000000,0.0146743226131101,0.0594271403285519;4.04700000000000,0.0146018927008991,0.0594271403285519;4.05700000000000,0.0145299977204853,0.0594271403285519;4.06700000000000,0.0144586324171492,0.0594271403285519;4.07700000000000,0.0143877916005356,0.0594271403285519;4.08700000000000,0.0143174701437100,0.0594271403285519;4.09700000000000,0.0142476629822305,0.0594271403285519;4.10700000000000,0.0141783651132367,0.0594271403285519;4.11700000000000,0.0141095715945525,0.0594271403285519;4.12700000000000,0.0140412775438056,0.0594271403285519;4.13700000000000,0.0139734781375606,0.0594271403285519;4.14700000000000,0.0139061686104679,0.0594271403285519;4.15700000000000,0.0138393442544258,0.0594271403285519;4.16700000000000,0.0137730004177578,0.0594271403285519;4.17700000000000,0.0137071325044022,0.0594271403285519;4.18700000000000,0.0136417359731170,0.0594271403285519;4.19700000000000,0.0135768063366966,0.0594271403285519;4.20700000000000,0.0135123391612019,0.0594271403285519;4.21700000000000,0.0134483300652041,0.0594271403285519;4.22700000000000,0.0133847747190392,0.0594271403285519;4.23700000000000,0.0133216688440767,0.0594271403285519;4.24700000000000,0.0132590082119990,0.0594271403285519;4.25700000000000,0.0131967886440931,0.0594271403285519;4.26700000000000,0.0131350060105543,0.0594271403285519;4.27700000000000,0.0130736562298003,0.0594271403285519;4.28700000000000,0.0130127352677977,0.0594271403285519;4.29700000000000,0.0129522391373982,0.0594271403285519;4.30700000000000,0.0128921638976867,0.0594271403285519;4.31700000000000,0.0128325056533395,0.0594271403285519;4.32700000000000,0.0127732605539924,0.0594271403285519;4.33700000000000,0.0127144247936196,0.0594271403285519;4.34700000000000,0.0126559946099229,0.0594271403285519;4.35700000000000,0.0125979662837295,0.0594271403285519;4.36700000000000,0.0125403361384006,0.0594271403285519;4.37700000000000,0.0124831005392489,0.0594271403285519;4.38700000000000,0.0124262558929653,0.0594271403285519;4.39700000000000,0.0123697986470555,0.0594271403285519;4.40700000000000,0.0123137252892841,0.0594271403285519;4.41700000000000,0.0122580323471292,0.0594271403285519;4.42700000000000,0.0122027163872440,0.0594271403285519;4.43700000000000,0.0121477740149283,0.0594271403285519;4.44700000000000,0.0120932018736073,0.0594271403285519;4.45700000000000,0.0120389966443188,0.0594271403285519;4.46700000000000,0.0119851550452090,0.0594271403285519;4.47700000000000,0.0119316738310354,0.0594271403285519;4.48700000000000,0.0118785497926779,0.0594271403285519;4.49700000000000,0.0118257797566574,0.0594271403285519;4.50700000000000,0.0117733605846616,0.0594271403285519;4.51700000000000,0.0117212891730789,0.0594271403285519;4.52700000000000,0.0116695624525385,0.0594271403285519;4.53700000000000,0.0116181773874584,0.0594271403285519;4.54700000000000,0.0115671309755998,0.0594271403285519;4.55700000000000,0.0115164202476287,0.0594271403285519;4.56700000000000,0.0114660422666841,0.0594271403285519;4.57700000000000,0.0114159941279528,0.0594271403285519;4.58700000000000,0.0113662729582503,0.0594271403285519;4.59700000000000,0.0113168759156090,0.0594271403285519;4.60700000000000,0.0112678001888718,0.0594271403285519;4.61700000000000,0.0112190429972919,0.0594271403285519;4.62700000000000,0.0111706015901394,0.0594271403285519;4.63700000000000,0.0111224732463129,0.0594271403285519;4.64700000000000,0.0110746552739575,0.0594271403285519;4.65700000000000,0.0110271450100884,0.0594271403285519;4.66700000000000,0.0109799398202204,0.0594271403285519;4.67700000000000,0.0109330370980019,0.0594271403285519;4.68700000000000,0.0108864342648562,0.0594271403285519;4.69700000000000,0.0108401287696263,0.0594271403285519;4.70700000000000,0.0107941180882258,0.0594271403285519;4.71700000000000,0.0107483997232955,0.0594271403285519;4.72700000000000,0.0107029712038638,0.0594271403285519;4.73700000000000,0.0106578300850134,0.0594271403285519;4.74700000000000,0.0106129739475519,0.0594271403285519;4.75700000000000,0.0105684003976879,0.0594271403285519;4.76700000000000,0.0105241070667117,0.0594271403285519;4.77700000000000,0.0104800916106804,0.0594271403285519;4.78700000000000,0.0104363517101077,0.0594271403285519;4.79700000000000,0.0103928850696589,0.0594271403285519;4.80700000000000,0.0103496894178489,0.0594271403285519;4.81700000000000,0.0103067625067463,0.0594271403285519;4.82700000000000,0.0102641021116802,0.0594271403285519;4.83700000000000,0.0102217060309526,0.0594271403285519;4.84700000000000,0.0101795720855538,0.0594271403285519;4.85700000000000,0.0101376981188827,0.0594271403285519;4.86700000000000,0.0100960819964708,0.0594271403285519;4.87700000000000,0.0100547216057100,0.0594271403285519;4.88700000000000,0.0100136148555848,0.0594271403285519;4.89700000000000,0.00997275967640769,0.0594271403285519;4.90700000000000,0.00993215401955896,0.0594271403285519;4.91700000000000,0.00989179585722958,0.0594271403285519;4.92700000000000,0.00985168318216818,0.0594271403285519;4.93700000000000,0.00981181400743143,0.0594271403285519;4.94700000000000,0.00977218636613800,0.0594271403285519;4.95700000000000,0.00973279831122597,0.0594271403285519;4.96700000000000,0.00969364791521361,0.0594271403285519;4.97700000000000,0.00965473326996365,0.0594271403285519;4.98700000000000,0.00961605248645071,0.0594271403285519;4.99700000000000,0.00957760369453209,0.0594271403285519], ...
        [0,0.156000000000000,0;0.149520685232551,0.363792000000000,0.00202099375946212;0.448562055697654,0.363792000000000,0.0181889438351591;0.458562055697654,0.355858678970057,0.0185944383170066;0.468562055697654,0.348263982074676,0.0189999327988542;0.478562055697654,0.340986681713556,0.0194054272807017;0.488562055697654,0.334007288251928,0.0198109217625492;0.498562055697654,0.327307875722739,0.0202164162443968;0.508562055697654,0.320871928092440,0.0206219107262443;0.518562055697654,0.314684203314530,0.0210274052080919;0.528562055697654,0.308730612815129,0.0214328996899394;0.538562055697654,0.302998114404798,0.0218383941717870;0.548562055697654,0.297474616903327,0.0222438886536345;0.558562055697654,0.292148895009601,0.0226493831354821;0.568562055697654,0.287010513155203,0.0230548776173296;0.578562055697654,0.282049757254799,0.0234603720991772;0.588562055697654,0.277257573414125,0.0238658665810247;0.598562055697654,0.272625512781900,0.0242713610628723;0.608562055697654,0.268145681838951,0.0246768555447198;0.618562055697654,0.263810697509262,0.0250823500265674;0.628562055697654,0.259613646555932,0.0254878445084149;0.638562055697654,0.255548048792340,0.0258933389902624;0.648562055697654,0.251607823696725,0.0262988334721100;0.658562055697654,0.247787260068440,0.0267043279539575;0.668562055697654,0.244080988407391,0.0271098224358051;0.678562055697654,0.240483955735760,0.0275153169176526;0.688562055697654,0.236991402613702,0.0279208113995002;0.698562055697654,0.233598842129193,0.0283263058813477;0.708562055697654,0.230302040666981,0.0287318003631953;0.718562055697654,0.227097000283331,0.0291372948450428;0.728562055697654,0.223979942532281,0.0295427893268904;0.738562055697654,0.220947293605838,0.0299482838087379;0.748562055697654,0.217995670665240,0.0303537782905855;0.758562055697654,0.215121869253373,0.0307592727724330;0.768562055697654,0.212322851689878,0.0311647672542806;0.778562055697654,0.209595736360585,0.0315702617361281;0.788562055697654,0.206937787821898,0.0319757562179756;0.798562055697654,0.204346407648680,0.0323812506998232;0.808562055697654,0.201819125961286,0.0327867451816708;0.818562055697654,0.199353593573649,0.0331922396635183;0.828562055697654,0.196947574709971,0.0335977341453659;0.838562055697654,0.194598940242530,0.0340032286272134;0.848562055697654,0.192305661407636,0.0344087231090609;0.858562055697654,0.190065803960741,0.0348142175909085;0.868562055697654,0.187877522735307,0.0352197120727560;0.878562055697654,0.185739056573277,0.0356252065546036;0.888562055697654,0.183648723597856,0.0360307010364511;0.898562055697654,0.181604916801949,0.0364361955182987;0.908562055697654,0.179606099927933,0.0368416900001462;0.918562055697654,0.177650803616553,0.0372471844819938;0.928562055697654,0.175737621804670,0.0376526789638413;0.938562055697654,0.173865208353286,0.0380581734456889;0.948562055697654,0.172032273888862,0.0384636679275364;0.958562055697654,0.170237582842348,0.0388691624093840;0.968562055697654,0.168479950671638,0.0392746568912315;0.978562055697654,0.166758241254328,0.0396801513730791;0.988562055697654,0.165071364438723,0.0400856458549266;0.998562055697654,0.163418273741987,0.0404911403367741;1.00856205569765,0.161797964185240,0.0408966348186217;1.01856205569765,0.160209470256174,0.0413021293004693;1.02856205569765,0.158651863990527,0.0417076237823168;1.03856205569765,0.157124253164384,0.0421131182641643;1.04856205569765,0.155625779589924,0.0425186127460119;1.05856205569765,0.154155617507765,0.0429241072278594;1.06856205569765,0.152712972069573,0.0433296017097070;1.07856205569765,0.151297077905089,0.0437350961915545;1.08856205569765,0.149907197768140,0.0441405906734021;1.09856205569765,0.148542621256593,0.0445460851552496;1.10856205569765,0.147202663601601,0.0449515796370971;1.11856205569765,0.145886664521784,0.0453570741189447;1.12856205569765,0.144593987138336,0.0457625686007923;1.13856205569765,0.143324016947298,0.0461680630826398;1.14856205569765,0.142076160845520,0.0465735575644874;1.15856205569765,0.140849846207070,0.0469790520463349;1.16856205569765,0.139644520007059,0.0473845465281824;1.17856205569765,0.138459647990079,0.0477900410100300;1.18856205569765,0.137294713880611,0.0481955354918775;1.19856205569765,0.136149218632969,0.0486010299737251;1.20856205569765,0.135022679718471,0.0490065244555726;1.21856205569765,0.133914630447716,0.0494120189374202;1.22856205569765,0.132824619325961,0.0498175134192677;1.23856205569765,0.131752209439715,0.0502230079011153;1.24856205569765,0.130696977872822,0.0506285023829628;1.25856205569765,0.129658515150375,0.0510339968648104;1.26856205569765,0.128636424708933,0.0514394913466579;1.27856205569765,0.127630322391602,0.0518449858285055;1.28856205569765,0.126639835966619,0.0522504803103530;1.29856205569765,0.125664604668192,0.0526559747922006;1.30856205569765,0.124704278758381,0.0530614692740481;1.31856205569765,0.123758519108924,0.0534669637558956;1.32856205569765,0.122826996801945,0.0538724582377432;1.33856205569765,0.121909392748557,0.0542779527195908;1.34856205569765,0.121005397324442,0.0546834472014383;1.35856205569765,0.120114710021518,0.0550889416832858;1.36856205569765,0.119237039114879,0.0554944361651334;1.37856205569765,0.118372101344236,0.0558999306469809;1.38856205569765,0.117519621609113,0.0563054251288285;1.39856205569765,0.116679332677132,0.0567109196106760;1.40856205569765,0.115850974904714,0.0571164140925236;1.41856205569765,0.115034295969595,0.0575219085743711;1.42856205569765,0.114229050614584,0.0579274030562187;1.43856205569765,0.113435000401997,0.0583328975380662;1.44856205569765,0.112651913478273,0.0587383920199137;1.45856205569765,0.111879564348263,0.0591438865017613;1.46856205569765,0.111117733658752,0.0595493809836088;1.47856205569765,0.110366207990752,0.0599548754654564;1.48856205569765,0.109624779660181,0.0603603699473039;1.49856205569765,0.108893246526512,0.0607658644291515;1.50856205569765,0.108171411809039,0.0611713589109990;1.51856205569765,0.107459083910398,0.0615768533928466;1.52856205569765,0.106756076247020,0.0619823478746941;1.53856205569765,0.106062207086191,0.0623878423565417;1.54856205569765,0.105377299389429,0.0627933368383892;1.55856205569765,0.104701180661886,0.0631988313202368;1.56856205569765,0.104033682807520,0.0636043258020843;1.57856205569765,0.103374641989758,0.0640098202839319;1.58856205569765,0.102723898497434,0.0644153147657794;1.59856205569765,0.102081296615753,0.0648208092476270;1.60856205569765,0.101446684502070,0.0652263037294745;1.61856205569765,0.100819914066269,0.0656317982113221;1.62856205569765,0.100200840855558,0.0660372926931696;1.63856205569765,0.0995893239434755,0.0664427871750171;1.64856205569765,0.0989852258229390,0.0668482816568647;1.65856205569765,0.0983884123031621,0.0672537761387122;1.66856205569765,0.0977987524102789,0.0676592706205598;1.67856205569765,0.0972161182915204,0.0680647651024073;1.68856205569765,0.0966403851227957,0.0684702595842549;1.69856205569765,0.0960714310195374,0.0688757540661024;1.70856205569765,0.0955091369506790,0.0692812485479500;1.71856205569765,0.0949533866556342,0.0696867430297975;1.72856205569765,0.0944040665641590,0.0700922375116451;1.73856205569765,0.0938610657189791,0.0704977319934926;1.74856205569765,0.0933242757010719,0.0709032264753402;1.75856205569765,0.0927935905574984,0.0713087209571877;1.76856205569765,0.0922689067316832,0.0717142154390352;1.77856205569765,0.0917501229960464,0.0721197099208828;1.78856205569765,0.0912371403868953,0.0725252044027303;1.79856205569765,0.0907298621414888,0.0729306988845779;1.80856205569765,0.0902281936371892,0.0733361933664255;1.81856205569765,0.0897320423326214,0.0737416878482730;1.82856205569765,0.0892413177107633,0.0741471823301205;1.83856205569765,0.0887559312238933,0.0745526768119681;1.84856205569765,0.0882757962403241,0.0749581712938156;1.85856205569765,0.0878008279928572,0.0753636657756632;1.86856205569765,0.0873309435288913,0.0757691602575107;1.87856205569765,0.0868660616621251,0.0761746547393583;1.88856205569765,0.0864061029257941,0.0765801492212058;1.89856205569765,0.0859509895273856,0.0769856437030533;1.90856205569765,0.0855006453047768,0.0773911381849009;1.91856205569765,0.0850549956837450,0.0777966326667485;1.92856205569765,0.0846139676367995,0.0782021271485960;1.93856205569765,0.0841774896432882,0.0786076216304436;1.94856205569765,0.0837454916507319,0.0790131161122911;1.95856205569765,0.0833179050373432,0.0794186105941386;1.96856205569765,0.0828946625756886,0.0798241050759862;1.97856205569765,0.0824756983974512,0.0802295995578337;1.98856205569765,0.0820609479592583,0.0806350940396813;1.99856205569765,0.0816503480095328,0.0810405885215288;2.00856205569765,0.0812438365563372,0.0814460830033764;2.01700000000000,0.0808638489867813,0.0817476875404655;2.02700000000000,0.0800679498040838,0.0817476875404655;2.03700000000000,0.0792837434633033,0.0817476875404655;2.04700000000000,0.0785110020360759,0.0817476875404655;2.05700000000000,0.0777495031208476,0.0817476875404655;2.06700000000000,0.0769990296828344,0.0817476875404655;2.07700000000000,0.0762593698993617,0.0817476875404655;2.08700000000000,0.0755303170103799,0.0817476875404655;2.09700000000000,0.0748116691739575,0.0817476875404655;2.10700000000000,0.0741032293265636,0.0817476875404655;2.11700000000000,0.0734048050479592,0.0817476875404655;2.12700000000000,0.0727162084305252,0.0817476875404655;2.13700000000000,0.0720372559528593,0.0817476875404655;2.14700000000000,0.0713677683574862,0.0817476875404655;2.15700000000000,0.0707075705325253,0.0817476875404655;2.16700000000000,0.0700564913971739,0.0817476875404655;2.17700000000000,0.0694143637908622,0.0817476875404655;2.18700000000000,0.0687810243659500,0.0817476875404655;2.19700000000000,0.0681563134838324,0.0817476875404655;2.20700000000000,0.0675400751143350,0.0817476875404655;2.21700000000000,0.0669321567382775,0.0817476875404655;2.22700000000000,0.0663324092530931,0.0817476875404655;2.23700000000000,0.0657406868813950,0.0817476875404655;2.24700000000000,0.0651568470823846,0.0817476875404655;2.25700000000000,0.0645807504660013,0.0817476875404655;2.26700000000000,0.0640122607097175,0.0817476875404655;2.27700000000000,0.0634512444778856,0.0817476875404655;2.28700000000000,0.0628975713435483,0.0817476875404655;2.29700000000000,0.0623511137126265,0.0817476875404655;2.30700000000000,0.0618117467504026,0.0817476875404655;2.31700000000000,0.0612793483102198,0.0817476875404655;2.32700000000000,0.0607537988643218,0.0817476875404655;2.33700000000000,0.0602349814367599,0.0817476875404655;2.34700000000000,0.0597227815382960,0.0817476875404655;2.35700000000000,0.0592170871032357,0.0817476875404655;2.36700000000000,0.0587177884281250,0.0817476875404655;2.37700000000000,0.0582247781122490,0.0817476875404655;2.38700000000000,0.0577379509998709,0.0817476875404655;2.39700000000000,0.0572572041241552,0.0817476875404655;2.40700000000000,0.0567824366527181,0.0817476875404655;2.41700000000000,0.0563135498347509,0.0817476875404655;2.42700000000000,0.0558504469496667,0.0817476875404655;2.43700000000000,0.0553930332572174,0.0817476875404655;2.44700000000000,0.0549412159490363,0.0817476875404655;2.45700000000000,0.0544949041015575,0.0817476875404655;2.46700000000000,0.0540540086302687,0.0817476875404655;2.47700000000000,0.0536184422452544,0.0817476875404655;2.48700000000000,0.0531881194079876,0.0817476875404655;2.49700000000000,0.0527629562893307,0.0817476875404655;2.50700000000000,0.0523428707287061,0.0817476875404655;2.51700000000000,0.0519277821944008,0.0817476875404655;2.52700000000000,0.0515176117449681,0.0817476875404655;2.53700000000000,0.0511122819916918,0.0817476875404655;2.54700000000000,0.0507117170620807,0.0817476875404655;2.55700000000000,0.0503158425643599,0.0817476875404655;2.56700000000000,0.0499245855529288,0.0817476875404655;2.57700000000000,0.0495378744947557,0.0817476875404655;2.58700000000000,0.0491556392366793,0.0817476875404655;2.59700000000000,0.0487778109735906,0.0817476875404655;2.60700000000000,0.0484043222174673,0.0817476875404655;2.61700000000000,0.0480351067672343,0.0817476875404655;2.62700000000000,0.0476700996794269,0.0817476875404655;2.63700000000000,0.0473092372396298,0.0817476875404655;2.64700000000000,0.0469524569346717,0.0817476875404655;2.65700000000000,0.0465996974255495,0.0817476875404655;2.66700000000000,0.0462508985210627,0.0817476875404655;2.67700000000000,0.0459060011521357,0.0817476875404655;2.68700000000000,0.0455649473468076,0.0817476875404655;2.69700000000000,0.0452276802058706,0.0817476875404655;2.70700000000000,0.0448941438791361,0.0817476875404655;2.71700000000000,0.0445642835423121,0.0817476875404655;2.72700000000000,0.0442380453744729,0.0817476875404655;2.73700000000000,0.0439153765361032,0.0817476875404655;2.74700000000000,0.0435962251477017,0.0817476875404655;2.75700000000000,0.0432805402689265,0.0817476875404655;2.76700000000000,0.0429682718782668,0.0817476875404655;2.77700000000000,0.0426593708532267,0.0817476875404655;2.78700000000000,0.0423537889510056,0.0817476875404655;2.79700000000000,0.0420514787896608,0.0817476875404655;2.80700000000000,0.0417523938297398,0.0817476875404655;2.81700000000000,0.0414564883563676,0.0817476875404655;2.82700000000000,0.0411637174617772,0.0817476875404655;2.83700000000000,0.0408740370282697,0.0817476875404655;2.84700000000000,0.0405874037115935,0.0817476875404655;2.85700000000000,0.0403037749247295,0.0817476875404655;2.86700000000000,0.0400231088220714,0.0817476875404655;2.87700000000000,0.0397453642839907,0.0817476875404655;2.88700000000000,0.0394705009017746,0.0817476875404655;2.89700000000000,0.0391984789629284,0.0817476875404655;2.90700000000000,0.0389292594368295,0.0817476875404655;2.91700000000000,0.0386628039607267,0.0817476875404655;2.92700000000000,0.0383990748260728,0.0817476875404655;2.93700000000000,0.0381380349651829,0.0817476875404655;2.94700000000000,0.0378796479382084,0.0817476875404655;2.95700000000000,0.0376238779204196,0.0817476875404655;2.96700000000000,0.0373706896897877,0.0817476875404655;2.97700000000000,0.0371200486148574,0.0817476875404655;2.98700000000000,0.0368719206429046,0.0817476875404655;2.99700000000000,0.0366262722883693,0.0817476875404655;3.00700000000000,0.0363830706215575,0.0817476875404655;3.01700000000000,0.0361422832576051,0.0817476875404655;3.02700000000000,0.0359038783456963,0.0817476875404655;3.03700000000000,0.0356678245585299,0.0817476875404655;3.04700000000000,0.0354340910820279,0.0817476875404655;3.05700000000000,0.0352026476052787,0.0817476875404655;3.06700000000000,0.0349734643107097,0.0817476875404655;3.07700000000000,0.0347465118644831,0.0817476875404655;3.08700000000000,0.0345217614071091,0.0817476875404655;3.09700000000000,0.0342991845442712,0.0817476875404655;3.10700000000000,0.0340787533378571,0.0817476875404655;3.11700000000000,0.0338604402971918,0.0817476875404655;3.12700000000000,0.0336442183704657,0.0817476875404655;3.13700000000000,0.0334300609363540,0.0817476875404655;3.14700000000000,0.0332179417958224,0.0817476875404655;3.15700000000000,0.0330078351641142,0.0817476875404655;3.16700000000000,0.0327997156629135,0.0817476875404655;3.17700000000000,0.0325935583126819,0.0817476875404655;3.18700000000000,0.0323893385251627,0.0817476875404655;3.19700000000000,0.0321870320960487,0.0817476875404655;3.20700000000000,0.0319866151978102,0.0817476875404655;3.21700000000000,0.0317880643726789,0.0817476875404655;3.22700000000000,0.0315913565257833,0.0817476875404655;3.23700000000000,0.0313964689184325,0.0817476875404655;3.24700000000000,0.0312033791615452,0.0817476875404655;3.25700000000000,0.0310120652092183,0.0817476875404655;3.26700000000000,0.0308225053524348,0.0817476875404655;3.27700000000000,0.0306346782129043,0.0817476875404655;3.28700000000000,0.0304485627370357,0.0817476875404655;3.29700000000000,0.0302641381900370,0.0817476875404655;3.30700000000000,0.0300813841501399,0.0817476875404655;3.31700000000000,0.0299002805029465,0.0817476875404655;3.32700000000000,0.0297208074358941,0.0817476875404655;3.33700000000000,0.0295429454328363,0.0817476875404655;3.34700000000000,0.0293666752687376,0.0817476875404655;3.35700000000000,0.0291919780044777,0.0817476875404655;3.36700000000000,0.0290188349817644,0.0817476875404655;3.37700000000000,0.0288472278181511,0.0817476875404655;3.38700000000000,0.0286771384021578,0.0817476875404655;3.39700000000000,0.0285085488884921,0.0817476875404655;3.40700000000000,0.0283414416933682,0.0817476875404655;3.41700000000000,0.0281757994899218,0.0817476875404655;3.42700000000000,0.0280116052037186,0.0817476875404655;3.43700000000000,0.0278488420083540,0.0817476875404655;3.44700000000000,0.0276874933211419,0.0817476875404655;3.45700000000000,0.0275275427988910,0.0817476875404655;3.46700000000000,0.0273689743337660,0.0817476875404655;3.47700000000000,0.0272117720492323,0.0817476875404655;3.48700000000000,0.0270559202960814,0.0817476875404655;3.49700000000000,0.0269014036485363,0.0817476875404655;3.50700000000000,0.0267482069004346,0.0817476875404655;3.51700000000000,0.0265963150614868,0.0817476875404655;3.52700000000000,0.0264457133536095,0.0817476875404655;3.53700000000000,0.0262963872073304,0.0817476875404655;3.54700000000000,0.0261483222582650,0.0817476875404655;3.55700000000000,0.0260015043436612,0.0817476875404655;3.56700000000000,0.0258559194990135,0.0817476875404655;3.57700000000000,0.0257115539547412,0.0817476875404655;3.58700000000000,0.0255683941329334,0.0817476875404655;3.59700000000000,0.0254264266441557,0.0817476875404655;3.60700000000000,0.0252856382843193,0.0817476875404655;3.61700000000000,0.0251460160316112,0.0817476875404655;3.62700000000000,0.0250075470434827,0.0817476875404655;3.63700000000000,0.0248702186536961,0.0817476875404655;3.64700000000000,0.0247340183694283,0.0817476875404655;3.65700000000000,0.0245989338684291,0.0817476875404655;3.66700000000000,0.0244649529962345,0.0817476875404655;3.67700000000000,0.0243320637634323,0.0817476875404655;3.68700000000000,0.0242002543429798,0.0817476875404655;3.69700000000000,0.0240695130675724,0.0817476875404655;3.70700000000000,0.0239398284270613,0.0817476875404655;3.71700000000000,0.0238111890659204,0.0817476875404655;3.72700000000000,0.0236835837807605,0.0817476875404655;3.73700000000000,0.0235570015178895,0.0817476875404655;3.74700000000000,0.0234314313709189,0.0817476875404655;3.75700000000000,0.0233068625784142,0.0817476875404655;3.76700000000000,0.0231832845215896,0.0817476875404655;3.77700000000000,0.0230606867220444,0.0817476875404655;3.78700000000000,0.0229390588395420,0.0817476875404655;3.79700000000000,0.0228183906698296,0.0817476875404655;3.80700000000000,0.0226986721424978,0.0817476875404655;3.81700000000000,0.0225798933188792,0.0817476875404655;3.82700000000000,0.0224620443899860,0.0817476875404655;3.83700000000000,0.0223451156744848,0.0817476875404655;3.84700000000000,0.0222290976167078,0.0817476875404655;3.85700000000000,0.0221139807847009,0.0817476875404655;3.86700000000000,0.0219997558683067,0.0817476875404655;3.87700000000000,0.0218864136772816,0.0817476875404655;3.88700000000000,0.0217739451394474,0.0817476875404655;3.89700000000000,0.0216623412988761,0.0817476875404655;3.90700000000000,0.0215515933141066,0.0817476875404655;3.91700000000000,0.0214416924563935,0.0817476875404655;3.92700000000000,0.0213326301079877,0.0817476875404655;3.93700000000000,0.0212243977604461,0.0817476875404655;3.94700000000000,0.0211169870129728,0.0817476875404655;3.95700000000000,0.0210103895707886,0.0817476875404655;3.96700000000000,0.0209045972435298,0.0817476875404655;3.97700000000000,0.0207996019436745,0.0817476875404655;3.98700000000000,0.0206953956849970,0.0817476875404655;3.99700000000000,0.0205919705810496,0.0817476875404655;4.00700000000000,0.0204893188436697,0.0817476875404655;4.01700000000000,0.0203874327815140,0.0817476875404655;4.02700000000000,0.0202863047986177,0.0817476875404655;4.03700000000000,0.0201859273929788,0.0817476875404655;4.04700000000000,0.0200862931551663,0.0817476875404655;4.05700000000000,0.0199873947669534,0.0817476875404655;4.06700000000000,0.0198892249999733,0.0817476875404655;4.07700000000000,0.0197917767143984,0.0817476875404655;4.08700000000000,0.0196950428576422,0.0817476875404655;4.09700000000000,0.0195990164630831,0.0817476875404655;4.10700000000000,0.0195036906488102,0.0817476875404655;4.11700000000000,0.0194090586163902,0.0817476875404655;4.12700000000000,0.0193151136496549,0.0817476875404655;4.13700000000000,0.0192218491135103,0.0817476875404655;4.14700000000000,0.0191292584527642,0.0817476875404655;4.15700000000000,0.0190373351909748,0.0817476875404655;4.16700000000000,0.0189460729293180,0.0817476875404655;4.17700000000000,0.0188554653454740,0.0817476875404655;4.18700000000000,0.0187655061925318,0.0817476875404655;4.19700000000000,0.0186761892979131,0.0817476875404655;4.20700000000000,0.0185875085623129,0.0817476875404655;4.21700000000000,0.0184994579586586,0.0817476875404655;4.22700000000000,0.0184120315310854,0.0817476875404655;4.23700000000000,0.0183252233939299,0.0817476875404655;4.24700000000000,0.0182390277307387,0.0817476875404655;4.25700000000000,0.0181534387932945,0.0817476875404655;4.26700000000000,0.0180684509006576,0.0817476875404655;4.27700000000000,0.0179840584382234,0.0817476875404655;4.28700000000000,0.0179002558567947,0.0817476875404655;4.29700000000000,0.0178170376716697,0.0817476875404655;4.30700000000000,0.0177343984617448,0.0817476875404654;4.31700000000000,0.0176523328686312,0.0817476875404655;4.32700000000000,0.0175708355957865,0.0817476875404655;4.33700000000000,0.0174899014076603,0.0817476875404655;4.34700000000000,0.0174095251288530,0.0817476875404655;4.35700000000000,0.0173297016432885,0.0817476875404655;4.36700000000000,0.0172504258934003,0.0817476875404655;4.37700000000000,0.0171716928793299,0.0817476875404655;4.38700000000000,0.0170934976581390,0.0817476875404655;4.39700000000000,0.0170158353430331,0.0817476875404655;4.40700000000000,0.0169387011025986,0.0817476875404655;4.41700000000000,0.0168620901600508,0.0817476875404655;4.42700000000000,0.0167859977924946,0.0817476875404655;4.43700000000000,0.0167104193301967,0.0817476875404655;4.44700000000000,0.0166353501558689,0.0817476875404655;4.45700000000000,0.0165607857039630,0.0817476875404655;4.46700000000000,0.0164867214599766,0.0817476875404655;4.47700000000000,0.0164131529597699,0.0817476875404654;4.48700000000000,0.0163400757888926,0.0817476875404655;4.49700000000000,0.0162674855819222,0.0817476875404655;4.50700000000000,0.0161953780218117,0.0817476875404655;4.51700000000000,0.0161237488392476,0.0817476875404655;4.52700000000000,0.0160525938120185,0.0817476875404655;4.53700000000000,0.0159819087643922,0.0817476875404655;4.54700000000000,0.0159116895665037,0.0817476875404655;4.55700000000000,0.0158419321337514,0.0817476875404655;4.56700000000000,0.0157726324262035,0.0817476875404655;4.57700000000000,0.0157037864480128,0.0817476875404655;4.58700000000000,0.0156353902468409,0.0817476875404655;4.59700000000000,0.0155674399132907,0.0817476875404655;4.60700000000000,0.0154999315803479,0.0817476875404655;4.61700000000000,0.0154328614228309,0.0817476875404655;4.62700000000000,0.0153662256568487,0.0817476875404655;4.63700000000000,0.0153000205392674,0.0817476875404655;4.64700000000000,0.0152342423671845,0.0817476875404655;4.65700000000000,0.0151688874774109,0.0817476875404655;4.66700000000000,0.0151039522459613,0.0817476875404655;4.67700000000000,0.0150394330875513,0.0817476875404655;4.68700000000000,0.0149753264551030,0.0817476875404655;4.69700000000000,0.0149116288392572,0.0817476875404655;4.70700000000000,0.0148483367678929,0.0817476875404655;4.71700000000000,0.0147854468056547,0.0817476875404655;4.72700000000000,0.0147229555534862,0.0817476875404655;4.73700000000000,0.0146608596481707,0.0817476875404655;4.74700000000000,0.0145991557618790,0.0817476875404655;4.75700000000000,0.0145378406017232,0.0817476875404655;4.76700000000000,0.0144769109093175,0.0817476875404655;4.77700000000000,0.0144163634603454,0.0817476875404655;4.78700000000000,0.0143561950641329,0.0817476875404655;4.79700000000000,0.0142964025632284,0.0817476875404655;4.80700000000000,0.0142369828329882,0.0817476875404655;4.81700000000000,0.0141779327811685,0.0817476875404655;4.82700000000000,0.0141192493475230,0.0817476875404655;4.83700000000000,0.0140609295034064,0.0817476875404655;4.84700000000000,0.0140029702513834,0.0817476875404655;4.85700000000000,0.0139453686248440,0.0817476875404655;4.86700000000000,0.0138881216876236,0.0817476875404655;4.87700000000000,0.0138312265336288,0.0817476875404655;4.88700000000000,0.0137746802864686,0.0817476875404655;4.89700000000000,0.0137184800990910,0.0817476875404655;4.90700000000000,0.0136626231534246,0.0817476875404655;4.91700000000000,0.0136071066600249,0.0817476875404655;4.92700000000000,0.0135519278577268,0.0817476875404655;4.93700000000000,0.0134970840133006,0.0817476875404655;4.94700000000000,0.0134425724211137,0.0817476875404655;4.95700000000000,0.0133883904027973,0.0817476875404655;4.96700000000000,0.0133345353069167,0.0817476875404655;4.97700000000000,0.0132810045086475,0.0817476875404655;4.98700000000000,0.0132277954094555,0.0817476875404655;4.99700000000000,0.0131749054367815,0.0817476875404655], ...
        [0,0.347077173600000,0;0.172222065078601,0.820490438390400,0.00604728420861178;0.516666195235803,0.820490438390400,0.0544255578775060;0.526666195235803,0.804911492070843,0.0554789567330798;0.536666195235803,0.789913128111712,0.0565323555886536;0.546666195235803,0.775463485258435,0.0575857544442274;0.556666195235803,0.761532991689844,0.0586391532998012;0.566666195235803,0.748094163009179,0.0596925521553750;0.576666195235803,0.735121421253382,0.0607459510109488;0.586666195235803,0.722590932412826,0.0617993498665226;0.596666195235803,0.710480460289845,0.0628527487220964;0.606666195235803,0.698769234810837,0.0639061475776701;0.616666195235803,0.687437833151247,0.0649595464332440;0.626666195235803,0.676468072242210,0.0660129452888178;0.636666195235803,0.665842911407470,0.0670663441443915;0.646666195235803,0.655546364033989,0.0681197429999653;0.656666195235803,0.645563417313264,0.0691731418555391;0.666666195235803,0.635879959205943,0.0702265407111129;0.676666195235803,0.626482711882478,0.0712799395666867;0.686666195235803,0.617359170979647,0.0723333384222605;0.696666195235803,0.608497550088588,0.0733867372778343;0.706666195235803,0.599886729956099,0.0744401361334081;0.716666195235803,0.591516211938869,0.0754935349889819;0.726666195235803,0.583376075300933,0.0765469338445557;0.736666195235803,0.575456937989166,0.0776003327001295;0.746666195235803,0.567749920560749,0.0786537315557033;0.756666195235803,0.560246612971016,0.0797071304112770;0.766666195235803,0.552939043960506,0.0807605292668508;0.776666195235803,0.545819652806980,0.0818139281224246;0.786666195235803,0.538881263231928,0.0828673269779984;0.796666195235803,0.532117059272296,0.0839207258335722;0.806666195235803,0.525520562946864,0.0849741246891460;0.816666195235803,0.519085613563473,0.0860275235447198;0.826666195235803,0.512806348528142,0.0870809224002936;0.836666195235803,0.506677185530423,0.0881343212558674;0.846666195235803,0.500692805991220,0.0891877201114412;0.856666195235803,0.494848139669895,0.0902411189670150;0.866666195235803,0.489138350337045,0.0912945178225888;0.876666195235803,0.483558822427845,0.0923479166781625;0.886666195235803,0.478105148598549,0.0934013155337364;0.896666195235803,0.472773118115647,0.0944547143893102;0.906666195235803,0.467558706013377,0.0955081132448839;0.916666195235803,0.462458062960940,0.0965615121004577;0.926666195235803,0.457467505785783,0.0976149109560315;0.936666195235803,0.452583508603941,0.0986683098116053;0.946666195235803,0.447802694512537,0.0997217086671791;0.956666195235803,0.443121827803307,0.100775107522753;0.966666195235803,0.438537806659429,0.101828506378327;0.976666195235803,0.434047656301009,0.102881905233900;0.986666195235803,0.429648522547397,0.103935304089474;0.996666195235803,0.425337665767051,0.104988702945048;1.00666619523580,0.421112455188012,0.106042101800622;1.01666619523580,0.416970363544153,0.107095500656196;1.02666619523580,0.412908962034304,0.108148899511769;1.03666619523580,0.408925915573140,0.109202298367343;1.04666619523580,0.405018978314303,0.110255697222917;1.05666619523580,0.401185989427743,0.111309096078491;1.06666619523580,0.397424869114569,0.112362494934065;1.07666619523580,0.393733614843996,0.113415893789638;1.08666619523580,0.390110297798060,0.114469292645212;1.09666619523580,0.386553059510851,0.115522691500786;1.10666619523580,0.383060108689953,0.116576090356360;1.11666619523580,0.379629718208677,0.117629489211934;1.12666619523580,0.376260222258466,0.118682888067507;1.13666619523580,0.372950013651617,0.119736286923081;1.14666619523580,0.369697541265135,0.120789685778655;1.15666619523580,0.366501307617192,0.121843084634229;1.16666619523580,0.363359866568211,0.122896483489803;1.17666619523580,0.360271821139190,0.123949882345376;1.18666619523580,0.357235821440322,0.125003281200950;1.19666619523580,0.354250562703487,0.126056680056524;1.20666619523580,0.351314783412560,0.127110078912098;1.21666619523580,0.348427263525937,0.128163477767671;1.22666619523580,0.345586822785993,0.129216876623245;1.23666619523580,0.342792319110568,0.130270275478819;1.24666619523580,0.340042647061863,0.131323674334393;1.25666619523580,0.337336736388440,0.132377073189967;1.26666619523580,0.334673550636288,0.133430472045540;1.27666619523580,0.332052085825164,0.134483870901114;1.28666619523580,0.329471369186655,0.135537269756688;1.29666619523580,0.326930457960642,0.136590668612262;1.30666619523580,0.324428438247017,0.137644067467836;1.31666619523580,0.321964423909739,0.138697466323409;1.32666619523580,0.319537555530444,0.139750865178983;1.33666619523580,0.317146999409033,0.140804264034557;1.34666619523580,0.314791946608785,0.141857662890131;1.35666619523580,0.312471612043700,0.142911061745705;1.36666619523580,0.310185233605915,0.143964460601278;1.37666619523580,0.307932071331143,0.145017859456852;1.38666619523580,0.305711406600228,0.146071258312426;1.39666619523580,0.303522541375001,0.147124657168000;1.40666619523580,0.301364797466723,0.148178056023574;1.41666619523580,0.299237515835523,0.149231454879147;1.42666619523580,0.297140055919288,0.150284853734721;1.43666619523580,0.295071794990586,0.151338252590295;1.44666619523580,0.293032127540262,0.152391651445869;1.45666619523580,0.291020464686421,0.153445050301443;1.46666619523580,0.289036233607585,0.154498449157016;1.47666619523580,0.287078876998894,0.155551848012590;1.48666619523580,0.285147852550239,0.156605246868164;1.49666619523580,0.283242632445329,0.157658645723738;1.50666619523580,0.281362702880699,0.158712044579311;1.51666619523580,0.279507563603747,0.159765443434885;1.52666619523580,0.277676727468932,0.160818842290459;1.53666619523580,0.275869720011296,0.161872241146033;1.54666619523580,0.274086079036526,0.162925640001607;1.55666619523580,0.272325354226831,0.163979038857181;1.56666619523580,0.270587106761896,0.165032437712754;1.57666619523580,0.268870908954272,0.166085836568328;1.58666619523580,0.267176343898549,0.167139235423902;1.59666619523580,0.265503005133717,0.168192634279476;1.60666619523580,0.263850496318127,0.169246033135049;1.61666619523580,0.262218430916527,0.170299431990623;1.62666619523580,0.260606431898631,0.171352830846197;1.63666619523580,0.259014131448745,0.172406229701771;1.64666619523580,0.257441170685974,0.173459628557345;1.65666619523580,0.255887199394556,0.174513027412918;1.66666619523580,0.254351875763909,0.175566426268492;1.67666619523580,0.252834866137982,0.176619825124066;1.68666619523580,0.251335844773517,0.177673223979640;1.69666619523580,0.249854493606863,0.178726622835214;1.70666619523580,0.248390502028988,0.179780021690787;1.71666619523580,0.246943566668355,0.180833420546361;1.72666619523580,0.245513391181340,0.181886819401935;1.73666619523580,0.244099686049894,0.182940218257509;1.74666619523580,0.242702168386155,0.183993617113082;1.75666619523580,0.241320561743729,0.185047015968656;1.76666619523580,0.239954595935392,0.186100414824230;1.77666619523580,0.238604006856933,0.187153813679804;1.78666619523580,0.237268536316923,0.188207212535378;1.79666619523580,0.235947931872168,0.189260611390951;1.80666619523580,0.234641946668623,0.190314010246525;1.81666619523580,0.233350339287565,0.191367409102099;1.82666619523580,0.232072873596810,0.192420807957673;1.83666619523580,0.230809318606802,0.193474206813247;1.84666619523580,0.229559448331372,0.194527605668820;1.85666619523580,0.228323041652991,0.195581004524394;1.86666619523580,0.227099882192367,0.196634403379968;1.87666619523580,0.225889758182199,0.197687802235542;1.88666619523580,0.224692462344957,0.198741201091116;1.89666619523580,0.223507791774514,0.199794599946689;1.90666619523580,0.222335547821519,0.200847998802263;1.91666619523580,0.221175535982347,0.201901397657837;1.92666619523580,0.220027565791510,0.202954796513411;1.93666619523580,0.218891450717406,0.204008195368985;1.94666619523580,0.217767008061274,0.205061594224558;1.95666619523580,0.216654058859251,0.206114993080132;1.96666619523580,0.215552427787419,0.207168391935706;1.97666619523580,0.214461943069732,0.208221790791280;1.98666619523580,0.213382436388720,0.209275189646854;1.99666619523580,0.212313742798886,0.210328588502427;2.00666619523580,0.211255700642682,0.211381987358001;2.01666619523580,0.210208151468992,0.212435386213575;2.02666619523580,0.209170939954027,0.213488785069149;2.03666619523580,0.208143913824544,0.214542183924723;2.04666619523580,0.207126923783330,0.215595582780296;2.05666619523580,0.206119823436841,0.216648981635870;2.06666619523580,0.205122469224961,0.217702380491444;2.07666619523580,0.204134720352776,0.218755779347018;2.08666619523580,0.203156438724316,0.219809178202591;2.09666619523580,0.202187488878194,0.220862577058165;2.10666619523580,0.201227737925074,0.221915975913739;2.11666619523580,0.200277055486918,0.222969374769313;2.12666619523580,0.199335313637936,0.224022773624887;2.13666619523580,0.198402386847207,0.225076172480460;2.14666619523580,0.197478151922897,0.226129571336034;2.15666619523580,0.196562487958028,0.227182970191608;2.16666619523580,0.195655276277751,0.228236369047182;2.17666619523580,0.194756400388071,0.229289767902756;2.18666619523580,0.193865745925985,0.230343166758329;2.19666619523580,0.192983200610969,0.231396565613903;2.20666619523580,0.192108654197797,0.232449964469477;2.21666619523580,0.191241998430634,0.233503363325051;2.22666619523580,0.190383126998356,0.234556762180625;2.23666619523580,0.189531935491086,0.235610161036198;2.24666619523580,0.188688321357873,0.236663559891772;2.25666619523580,0.187852183865513,0.237716958747346;2.26666619523580,0.187023424058443,0.238770357602920;2.27666619523580,0.186201944719708,0.239823756458493;2.28666619523580,0.185387650332938,0.240877155314067;2.29666619523580,0.184580447045331,0.241930554169641;2.30666619523580,0.183780242631591,0.242983953025215;2.31666619523580,0.182986946458799,0.244037351880789;2.32666619523580,0.182200469452199,0.245090750736362;2.33666619523580,0.181420724061848,0.246144149591936;2.34666619523580,0.180647624230138,0.247197548447510;2.35666619523580,0.179881085360121,0.248250947303084;2.36666619523580,0.179121024284663,0.249304346158658;2.37666619523580,0.178367359236354,0.250357745014232;2.38666619523580,0.177620009818189,0.251411143869805;2.39666619523580,0.176878896974977,0.252464542725379;2.40666619523580,0.176143942965463,0.253517941580953;2.41666619523580,0.175415071335146,0.254571340436527;2.42666619523580,0.174692206889762,0.255624739292100;2.43666619523580,0.173975275669428,0.256678138147674;2.44666619523580,0.173264204923413,0.257731537003248;2.45666619523580,0.172558923085533,0.258784935858822;2.46666619523580,0.171859359750134,0.259838334714396;2.47666619523580,0.171165445648667,0.260891733569969;2.48666619523580,0.170477112626822,0.261945132425543;2.49666619523580,0.169794293622214,0.262998531281117;2.50666619523580,0.169116922642604,0.264051930136691;2.51666619523580,0.168444934744635,0.265105328992265;2.52666619523580,0.167778266013078,0.266158727847838;2.53666619523580,0.167116853540565,0.267212126703412;2.54666619523580,0.166460635407803,0.268265525558986;2.55666619523580,0.165809550664249,0.269318924414560;2.56666619523580,0.165163539309239,0.270372323270134;2.57666619523580,0.164522542273556,0.271425722125707;2.58666619523580,0.163886501401422,0.272479120981281;2.59666619523580,0.163255359432916,0.273532519836855;2.60666619523580,0.162629059986783,0.274585918692429;2.61666619523580,0.162007547543649,0.275639317548002;2.62666619523580,0.161390767429611,0.276692716403576;2.63666619523580,0.160778665800208,0.277746115259150;2.64500000000000,0.160211493740046,0.278518657413711;2.65500000000000,0.159006900314966,0.278518657413711;2.66500000000000,0.157815841583828,0.278518657413711;2.67500000000000,0.156638115537456,0.278518657413711;2.68500000000000,0.155473523921462,0.278518657413711;2.69500000000000,0.154321872152816,0.278518657413711;2.70500000000000,0.153182969238551,0.278518657413711;2.71500000000000,0.152056627696578,0.278518657413711;2.72500000000000,0.150942663478523,0.278518657413711;2.73500000000000,0.149840895894536,0.278518657413711;2.74500000000000,0.148751147540015,0.278518657413711;2.75500000000000,0.147673244224190,0.278518657413711;2.76500000000000,0.146607014900504,0.278518657413711;2.77500000000000,0.145552291598760,0.278518657413711;2.78500000000000,0.144508909358961,0.278518657413711;2.79500000000000,0.143476706166801,0.278518657413711;2.80500000000000,0.142455522890777,0.278518657413711;2.81500000000000,0.141445203220846,0.278518657413711;2.82500000000000,0.140445593608609,0.278518657413711;2.83500000000000,0.139456543208969,0.278518657413711;2.84500000000000,0.138477903823216,0.278518657413711;2.85500000000000,0.137509529843511,0.278518657413711;2.86500000000000,0.136551278198722,0.278518657413711;2.87500000000000,0.135603008301575,0.278518657413711;2.88500000000000,0.134664581997087,0.278518657413711;2.89500000000000,0.133735863512244,0.278518657413711;2.90500000000000,0.132816719406887,0.278518657413711;2.91500000000000,0.131907018525778,0.278518657413711;2.92500000000000,0.131006631951810,0.278518657413711;2.93500000000000,0.130115432960331,0.278518657413711;2.94500000000000,0.129233296974551,0.278518657413711;2.95500000000000,0.128360101522007,0.278518657413711;2.96500000000000,0.127495726192050,0.278518657413711;2.97500000000000,0.126640052594331,0.278518657413711;2.98500000000000,0.125792964318264,0.278518657413711;2.99500000000000,0.124954346893426,0.278518657413711;3.00500000000000,0.124124087750887,0.278518657413711;3.01500000000000,0.123302076185430,0.278518657413711;3.02500000000000,0.122488203318648,0.278518657413711;3.03500000000000,0.121682362062886,0.278518657413711;3.04500000000000,0.120884447086015,0.278518657413711;3.05500000000000,0.120094354777010,0.278518657413711;3.06500000000000,0.119311983212314,0.278518657413711;3.07500000000000,0.118537232122964,0.278518657413711;3.08500000000000,0.117770002862463,0.278518657413711;3.09500000000000,0.117010198375378,0.278518657413711;3.10500000000000,0.116257723166645,0.278518657413711;3.11500000000000,0.115512483271562,0.278518657413711;3.12500000000000,0.114774386226453,0.278518657413711;3.13500000000000,0.114043341039985,0.278518657413711;3.14500000000000,0.113319258165125,0.278518657413711;3.15500000000000,0.112602049471717,0.278518657413711;3.16500000000000,0.111891628219662,0.278518657413711;3.17500000000000,0.111187909032694,0.278518657413711;3.18500000000000,0.110490807872726,0.278518657413711;3.19500000000000,0.109800242014758,0.278518657413711;3.20500000000000,0.109116130022338,0.278518657413711;3.21500000000000,0.108438391723546,0.278518657413711;3.22500000000000,0.107766948187508,0.278518657413711;3.23500000000000,0.107101721701416,0.278518657413711;3.24500000000000,0.106442635748035,0.278518657413711;3.25500000000000,0.105789614983703,0.278518657413711;3.26500000000000,0.105142585216795,0.278518657413711;3.27500000000000,0.104501473386652,0.278518657413711;3.28500000000000,0.103866207542953,0.278518657413711;3.29500000000000,0.103236716825530,0.278518657413711;3.30500000000000,0.102612931444605,0.278518657413711;3.31500000000000,0.101994782661444,0.278518657413711;3.32500000000000,0.101382202769423,0.278518657413711;3.33500000000000,0.100775125075487,0.278518657413711;3.34500000000000,0.100173483881992,0.278518657413711;3.35500000000000,0.0995772144689360,0.278518657413711;3.36500000000000,0.0989862530765489,0.278518657413711;3.37500000000000,0.0984005368882484,0.278518657413711;3.38500000000000,0.0978200040139468,0.278518657413711;3.39500000000000,0.0972445934737002,0.278518657413711;3.40500000000000,0.0966742451816953,0.278518657413711;3.41500000000000,0.0961088999305626,0.278518657413711;3.42500000000000,0.0955484993760098,0.278518657413711;3.43500000000000,0.0949929860217688,0.278518657413711;3.44500000000000,0.0944423032048470,0.278518657413711;3.45500000000000,0.0938963950810780,0.278518657413711;3.46500000000000,0.0933552066109626,0.278518657413711;3.47500000000000,0.0928186835457962,0.278518657413711;3.48500000000000,0.0922867724140726,0.278518657413711;3.49500000000000,0.0917594205081615,0.278518657413711;3.50500000000000,0.0912365758712502,0.278518657413711;3.51500000000000,0.0907181872845460,0.278518657413711;3.52500000000000,0.0902042042547320,0.278518657413711;3.53500000000000,0.0896945770016709,0.278518657413711;3.54500000000000,0.0891892564463511,0.278518657413711;3.55500000000000,0.0886881941990702,0.278518657413711;3.56500000000000,0.0881913425478505,0.278518657413711;3.57500000000000,0.0876986544470795,0.278518657413711;3.58500000000000,0.0872100835063738,0.278518657413711;3.59500000000000,0.0867255839796584,0.278518657413711;3.60500000000000,0.0862451107544580,0.278518657413711;3.61500000000000,0.0857686193413952,0.278518657413711;3.62500000000000,0.0852960658638919,0.278518657413711;3.63500000000000,0.0848274070480678,0.278518657413711;3.64500000000000,0.0843626002128330,0.278518657413711;3.65500000000000,0.0839016032601709,0.278518657413711;3.66500000000000,0.0834443746656049,0.278518657413711;3.67500000000000,0.0829908734688476,0.278518657413711;3.68500000000000,0.0825410592646270,0.278518657413711;3.69500000000000,0.0820948921936864,0.278518657413711;3.70500000000000,0.0816523329339536,0.278518657413711;3.71500000000000,0.0812133426918773,0.278518657413711;3.72500000000000,0.0807778831939249,0.278518657413711;3.73500000000000,0.0803459166782403,0.278518657413711;3.74500000000000,0.0799174058864569,0.278518657413711;3.75500000000000,0.0794923140556634,0.278518657413711;3.76500000000000,0.0790706049105185,0.278518657413711;3.77500000000000,0.0786522426555119,0.278518657413711;3.78500000000000,0.0782371919673679,0.278518657413711;3.79500000000000,0.0778254179875889,0.278518657413711;3.80500000000000,0.0774168863151366,0.278518657413711;3.81500000000000,0.0770115629992463,0.278518657413711;3.82500000000000,0.0766094145323733,0.278518657413711;3.83500000000000,0.0762104078432678,0.278518657413711;3.84500000000000,0.0758145102901750,0.278518657413711;3.85500000000000,0.0754216896541595,0.278518657413711;3.86500000000000,0.0750319141325495,0.278518657413711;3.87500000000000,0.0746451523325006,0.278518657413711;3.88500000000000,0.0742613732646737,0.278518657413711;3.89500000000000,0.0738805463370277,0.278518657413711;3.90500000000000,0.0735026413487226,0.278518657413711;3.91500000000000,0.0731276284841323,0.278518657413711;3.92500000000000,0.0727554783069629,0.278518657413711;3.93500000000000,0.0723861617544762,0.278518657413711;3.94500000000000,0.0720196501318159,0.278518657413711;3.95500000000000,0.0716559151064331,0.278518657413711;3.96500000000000,0.0712949287026110,0.278518657413711;3.97500000000000,0.0709366632960851,0.278518657413711;3.98500000000000,0.0705810916087590,0.278518657413711;3.99500000000000,0.0702281867035111,0.278518657413711;4.00500000000000,0.0698779219790932,0.278518657413711;4.01500000000000,0.0695302711651174,0.278518657413711;4.02500000000000,0.0691852083171300,0.278518657413711;4.03500000000000,0.0688427078117712,0.278518657413711;4.04500000000000,0.0685027443420179,0.278518657413711;4.05500000000000,0.0681652929125088,0.278518657413711;4.06500000000000,0.0678303288349502,0.278518657413711;4.07500000000000,0.0674978277236000,0.278518657413711;4.08500000000000,0.0671677654908294,0.278518657413711;4.09500000000000,0.0668401183427602,0.278518657413711;4.10500000000000,0.0665148627749769,0.278518657413711;4.11500000000000,0.0661919755683105,0.278518657413711;4.12500000000000,0.0658714337846952,0.278518657413711;4.13500000000000,0.0655532147630941,0.278518657413711;4.14500000000000,0.0652372961154940,0.278518657413711;4.15500000000000,0.0649236557229675,0.278518657413711;4.16500000000000,0.0646122717318018,0.278518657413711;4.17500000000000,0.0643031225496908,0.278518657413711;4.18500000000000,0.0639961868419930,0.278518657413711;4.19500000000000,0.0636914435280496,0.278518657413711;4.20500000000000,0.0633888717775653,0.278518657413711;4.21500000000000,0.0630884510070488,0.278518657413711;4.22500000000000,0.0627901608763113,0.278518657413711;4.23500000000000,0.0624939812850246,0.278518657413711;4.24500000000000,0.0621998923693338,0.278518657413711;4.25500000000000,0.0619078744985276,0.278518657413711;4.26500000000000,0.0616179082717616,0.278518657413711;4.27500000000000,0.0613299745148363,0.278518657413711;4.28500000000000,0.0610440542770270,0.278518657413711;4.29500000000000,0.0607601288279658,0.278518657413711;4.30500000000000,0.0604781796545737,0.278518657413711;4.31500000000000,0.0601981884580429,0.278518657413711;4.32500000000000,0.0599201371508680,0.278518657413711;4.33500000000000,0.0596440078539239,0.278518657413711;4.34500000000000,0.0593697828935925,0.278518657413711;4.35500000000000,0.0590974447989341,0.278518657413711;4.36500000000000,0.0588269762989050,0.278518657413711;4.37500000000000,0.0585583603196189,0.278518657413711;4.38500000000000,0.0582915799816522,0.278518657413711;4.39500000000000,0.0580266185973928,0.278518657413711;4.40500000000000,0.0577634596684299,0.278518657413711;4.41500000000000,0.0575020868829856,0.278518657413711;4.42500000000000,0.0572424841133878,0.278518657413711;4.43500000000000,0.0569846354135816,0.278518657413711;4.44500000000000,0.0567285250166808,0.278518657413711;4.45500000000000,0.0564741373325576,0.278518657413711;4.46500000000000,0.0562214569454701,0.278518657413711;4.47500000000000,0.0559704686117265,0.278518657413711;4.48500000000000,0.0557211572573861,0.278518657413711;4.49500000000000,0.0554735079759963,0.278518657413711;4.50500000000000,0.0552275060263638,0.278518657413711;4.51500000000000,0.0549831368303614,0.278518657413711;4.52500000000000,0.0547403859707679,0.278518657413711;4.53500000000000,0.0544992391891416,0.278518657413711;4.54500000000000,0.0542596823837268,0.278518657413711;4.55500000000000,0.0540217016073918,0.278518657413711;4.56500000000000,0.0537852830655989,0.278518657413711;4.57500000000000,0.0535504131144056,0.278518657413711;4.58500000000000,0.0533170782584957,0.278518657413711;4.59500000000000,0.0530852651492411,0.278518657413711;4.60500000000000,0.0528549605827921,0.278518657413711;4.61500000000000,0.0526261514981978,0.278518657413711;4.62500000000000,0.0523988249755538,0.278518657413711;4.63500000000000,0.0521729682341783,0.278518657413711;4.64500000000000,0.0519485686308161,0.278518657413711;4.65500000000000,0.0517256136578690,0.278518657413711;4.66500000000000,0.0515040909416526,0.278518657413711;4.67500000000000,0.0512839882406797,0.278518657413711;4.68500000000000,0.0510652934439692,0.278518657413711;4.69500000000000,0.0508479945693799,0.278518657413711;4.70500000000000,0.0506320797619691,0.278518657413711;4.71500000000000,0.0504175372923761,0.278518657413711;4.72500000000000,0.0502043555552288,0.278518657413711;4.73500000000000,0.0499925230675743,0.278518657413711;4.74500000000000,0.0497820284673326,0.278518657413711;4.75500000000000,0.0495728605117732,0.278518657413711;4.76500000000000,0.0493650080760136,0.278518657413711;4.77500000000000,0.0491584601515399,0.278518657413711;4.78500000000000,0.0489532058447497,0.278518657413711;4.79500000000000,0.0487492343755152,0.278518657413711;4.80500000000000,0.0485465350757678,0.278518657413711;4.81500000000000,0.0483450973881035,0.278518657413711;4.82500000000000,0.0481449108644078,0.278518657413711;4.83500000000000,0.0479459651645011,0.278518657413711;4.84500000000000,0.0477482500548033,0.278518657413711;4.85500000000000,0.0475517554070179,0.278518657413711;4.86500000000000,0.0473564711968348,0.278518657413711;4.87500000000000,0.0471623875026516,0.278518657413711;4.88500000000000,0.0469694945043139,0.278518657413711;4.89500000000000,0.0467777824818723,0.278518657413711;4.90500000000000,0.0465872418143588,0.278518657413711;4.91500000000000,0.0463978629785791,0.278518657413711;4.92500000000000,0.0462096365479226,0.278518657413711;4.93500000000000,0.0460225531911898,0.278518657413711;4.94500000000000,0.0458366036714355,0.278518657413711;4.95500000000000,0.0456517788448287,0.278518657413711;4.96500000000000,0.0454680696595283,0.278518657413711;4.97500000000000,0.0452854671545751,0.278518657413711;4.98500000000000,0.0451039624587989,0.278518657413711;4.99500000000000,0.0449235467897409,0.278518657413711], ...
        [0,0.407159360000000,0;0.177830386091840,0.977182464000000,0.00767886227286544;0.533491158275520,0.977182464000000,0.0691097604557889;0.543491158275520,0.959202733343469,0.0704051851199799;0.553491158275520,0.941872687162930,0.0717006097841708;0.563491158275520,0.925157736567336,0.0729960344483617;0.573491158275520,0.909025705179977,0.0742914591125526;0.583491158275520,0.893446622407471,0.0755868837767435;0.593491158275521,0.878392537608575,0.0768823084409344;0.603491158275521,0.863837352738616,0.0781777331051254;0.613491158275520,0.849756671361450,0.0794731577693163;0.623491158275520,0.836127662191349,0.0807685824335072;0.633491158275520,0.822928935559276,0.0820640070976981;0.643491158275520,0.810140431397624,0.0833594317618891;0.653491158275520,0.797743317509573,0.0846548564260800;0.663491158275520,0.785719897038033,0.0859502810902709;0.673491158275520,0.774053524177996,0.0872457057544618;0.683491158275520,0.762728527288044,0.0885411304186527;0.693491158275520,0.751730138654154,0.0898365550828436;0.703491158275520,0.741044430243875,0.0911319797470345;0.713491158275520,0.730658254863161,0.0924274044112254;0.723491158275520,0.720559192193138,0.0937228290754164;0.733491158275520,0.710735499241103,0.0950182537396073;0.743491158275520,0.701176064790133,0.0963136784037982;0.753491158275520,0.691870367475848,0.0976091030679891;0.763491158275520,0.682808437157774,0.0989045277321800;0.773491158275520,0.673980819287131,0.100199952396371;0.783491158275520,0.665378542003357,0.101495377060562;0.793491158275520,0.656993085718634,0.102790801724753;0.803491158275520,0.648816354973660,0.104086226388944;0.813491158275520,0.640840652369232,0.105381651053135;0.823491158275520,0.633058654397193,0.106677075717326;0.833491158275520,0.625463389011211,0.107972500381516;0.843491158275521,0.618048214793026,0.109267925045707;0.853491158275521,0.610806801583289,0.110563349709898;0.863491158275521,0.603733112458282,0.111858774374089;0.873491158275521,0.596821386944652,0.113154199038280;0.883491158275521,0.590066125374072,0.114449623702471;0.893491158275520,0.583462074288520,0.115745048366662;0.903491158275520,0.577004212814787,0.117040473030853;0.913491158275520,0.570687739933933,0.118335897695044;0.923491158275520,0.564508062577848,0.119631322359235;0.933491158275520,0.558460784490922,0.120926747023426;0.943491158275520,0.552541695800026,0.122222171687617;0.953491158275520,0.546746763240826,0.123517596351807;0.963491158275520,0.541072120992739,0.124813021015998;0.973491158275520,0.535514062078766,0.126108445680189;0.983491158275520,0.530069030289993,0.127403870344380;0.993491158275520,0.524733612597801,0.128699295008571;1.00349115827552,0.519504532019756,0.129994719672762;1.01349115827552,0.514378640907852,0.131290144336953;1.02349115827552,0.509352914630211,0.132585569001144;1.03349115827552,0.504424445619599,0.133880993665335;1.04349115827552,0.499590437764150,0.135176418329526;1.05349115827552,0.494848201117551,0.136471842993717;1.06349115827552,0.490195146907679,0.137767267657908;1.07349115827552,0.485628782824205,0.139062692322098;1.08349115827552,0.481146708567160,0.140358116986289;1.09349115827552,0.476746611639756,0.141653541650480;1.10349115827552,0.472426263369954,0.142948966314671;1.11349115827552,0.468183515146416,0.144244390978862;1.12349115827552,0.464016294855470,0.145539815643053;1.13349115827552,0.459922603506686,0.146835240307244;1.14349115827552,0.455900512035509,0.148130664971435;1.15349115827552,0.451948158272199,0.149426089635626;1.16349115827552,0.448063744067092,0.150721514299817;1.17349115827552,0.444245532562835,0.152016938964008;1.18349115827552,0.440491845604919,0.153312363628199;1.19349115827552,0.436801061282382,0.154607788292390;1.20349115827552,0.433171611591133,0.155903212956580;1.21349115827552,0.429601980212799,0.157198637620771;1.22349115827552,0.426090700402504,0.158494062284962;1.23349115827552,0.422636352979388,0.159789486949153;1.24349115827552,0.419237564414092,0.161084911613344;1.25349115827552,0.415893005007778,0.162380336277535;1.26349115827552,0.412601387157635,0.163675760941726;1.27349115827552,0.409361463704092,0.164971185605917;1.28349115827552,0.406172026355306,0.166266610270108;1.29349115827552,0.403031904184724,0.167562034934299;1.30349115827552,0.399939962197807,0.168857459598490;1.31349115827552,0.396895099964224,0.170152884262680;1.32349115827552,0.393896250312056,0.171448308926871;1.33349115827552,0.390942378080751,0.172743733591062;1.34349115827552,0.388032478929777,0.174039158255253;1.35349115827552,0.385165578200080,0.175334582919444;1.36349115827552,0.382340729825653,0.176630007583635;1.37349115827552,0.379557015292640,0.177925432247826;1.38349115827552,0.376813542643593,0.179220856912017;1.39349115827552,0.374109445524600,0.180516281576208;1.40349115827552,0.371443882273141,0.181811706240399;1.41349115827552,0.368816035044678,0.183107130904590;1.42349115827552,0.366225108976044,0.184402555568781;1.43349115827552,0.363670331383857,0.185697980232971;1.44349115827552,0.361150950996253,0.186993404897162;1.45349115827552,0.358666237216331,0.188288829561353;1.46349115827552,0.356215479415792,0.189584254225544;1.47349115827552,0.353797986257348,0.190879678889735;1.48349115827552,0.351413085044532,0.192175103553926;1.49349115827552,0.349060121097626,0.193470528218117;1.50349115827552,0.346738457154501,0.194765952882308;1.51349115827552,0.344447472795203,0.196061377546499;1.52349115827552,0.342186563889206,0.197356802210690;1.53349115827552,0.339955142064290,0.198652226874881;1.54349115827552,0.337752634196061,0.199947651539072;1.55349115827552,0.335578481917197,0.201243076203262;1.56349115827552,0.333432141145514,0.202538500867453;1.57349115827552,0.331313081630042,0.203833925531644;1.58349115827552,0.329220786514288,0.205129350195835;1.59349115827552,0.327154751915950,0.206424774860026;1.60349115827552,0.325114486522358,0.207720199524217;1.61349115827552,0.323099511200957,0.209015624188408;1.62349115827552,0.321109358624185,0.210311048852599;1.63349115827552,0.319143572908129,0.211606473516790;1.64349115827552,0.317201709264378,0.212901898180981;1.65349115827552,0.315283333664503,0.214197322845172;1.66349115827552,0.313388022516644,0.215492747509363;1.67349115827552,0.311515362353685,0.216788172173554;1.68349115827552,0.309664949532552,0.218083596837744;1.69349115827552,0.307836389944158,0.219379021501935;1.70349115827552,0.306029298733566,0.220674446166126;1.71349115827552,0.304243300029951,0.221969870830317;1.72349115827552,0.302478026685965,0.223265295494508;1.73349115827552,0.300733120026118,0.224560720158699;1.74349115827552,0.299008229603826,0.225856144822890;1.75349115827552,0.297303012966761,0.227151569487081;1.76349115827552,0.295617135430196,0.228446994151272;1.77349115827552,0.293950269858012,0.229742418815463;1.78349115827552,0.292302096451073,0.231037843479654;1.79349115827552,0.290672302542682,0.232333268143844;1.80349115827552,0.289060582400840,0.233628692808035;1.81349115827552,0.287466637037049,0.234924117472226;1.82349115827552,0.285890174021408,0.236219542136417;1.83349115827552,0.284330907303753,0.237514966800608;1.84349115827552,0.282788557040626,0.238810391464799;1.85349115827552,0.281262849427844,0.240105816128990;1.86349115827552,0.279753516538450,0.241401240793181;1.87349115827552,0.278260296165871,0.242696665457372;1.88349115827552,0.276782931672052,0.243992090121563;1.89349115827552,0.275321171840418,0.245287514785754;1.90349115827552,0.273874770733466,0.246582939449945;1.91349115827552,0.272443487554816,0.247878364114135;1.92349115827552,0.271027086515577,0.249173788778326;1.93349115827552,0.269625336704850,0.250469213442517;1.94349115827552,0.268238011964231,0.251764638106708;1.95349115827552,0.266864890766176,0.253060062770899;1.96349115827552,0.265505756096068,0.254355487435090;1.97349115827552,0.264160395337888,0.255650912099281;1.98349115827552,0.262828600163325,0.256946336763472;1.99349115827552,0.261510166424242,0.258241761427663;2.00349115827552,0.260204894048350,0.259537186091854;2.01349115827552,0.258912586937991,0.260832610756045;2.02349115827552,0.257633052871933,0.262128035420236;2.03349115827552,0.256366103410051,0.263423460084426;2.04349115827552,0.255111553800811,0.264718884748617;2.05349115827552,0.253869222891458,0.266014309412808;2.06349115827552,0.252638933040817,0.267309734076999;2.07349115827552,0.251420510034611,0.268605158741190;2.08349115827552,0.250213783003224,0.269900583405381;2.09349115827552,0.249018584341820,0.271196008069572;2.10349115827552,0.247834749632736,0.272491432733763;2.11349115827552,0.246662117570083,0.273786857397954;2.12349115827552,0.245500529886476,0.275082282062145;2.13349115827552,0.244349831281824,0.276377706726336;2.14349115827552,0.243209869354114,0.277673131390527;2.15349115827552,0.242080494532121,0.278968556054718;2.16349115827552,0.240961560009989,0.280263980718908;2.17349115827552,0.239852921683615,0.281559405383099;2.18349115827552,0.238754438088778,0.282854830047290;2.19349115827552,0.237665970340969,0.284150254711481;2.20349115827552,0.236587382076848,0.285445679375672;2.21349115827552,0.235518539397299,0.286741104039863;2.22349115827552,0.234459310812015,0.288036528704054;2.23349115827552,0.233409567185570,0.289331953368245;2.24349115827552,0.232369181684943,0.290627378032436;2.25349115827552,0.231338029728426,0.291922802696627;2.26349115827552,0.230315988935898,0.293218227360818;2.27349115827552,0.229302939080403,0.294513652025009;2.28349115827552,0.228298762041007,0.295809076689199;2.29349115827552,0.227303341756882,0.297104501353390;2.30349115827552,0.226316564182589,0.298399926017581;2.31349115827552,0.225338317244522,0.299695350681772;2.32349115827552,0.224368490798478,0.300990775345963;2.33349115827552,0.223406976588310,0.302286200010154;2.34349115827552,0.222453668205646,0.303581624674345;2.35349115827552,0.221508461050635,0.304877049338536;2.36349115827552,0.220571252293686,0.306172474002727;2.37349115827552,0.219641940838177,0.307467898666918;2.38349115827552,0.218720427284096,0.308763323331109;2.39349115827552,0.217806613892608,0.310058747995299;2.40349115827552,0.216900404551489,0.311354172659490;2.41349115827552,0.216001704741433,0.312649597323681;2.42349115827552,0.215110421503184,0.313945021987872;2.43349115827552,0.214226463405487,0.315240446652063;2.44349115827552,0.213349740513816,0.316535871316254;2.45349115827552,0.212480164359877,0.317831295980445;2.46349115827552,0.211617647911843,0.319126720644636;2.47349115827552,0.210762105545322,0.320422145308827;2.48349115827552,0.209913453015020,0.321717569973018;2.49349115827552,0.209071607427085,0.323012994637209;2.50349115827552,0.208236487212116,0.324308419301400;2.51349115827552,0.207408012098820,0.325603843965591;2.52349115827552,0.206586103088287,0.326899268629781;2.53349115827552,0.205770682428880,0.328194693293972;2.54349115827552,0.204961673591718,0.329490117958163;2.55349115827552,0.204159001246730,0.330785542622354;2.56349115827552,0.203362591239278,0.332080967286545;2.57349115827552,0.202572370567311,0.333376391950736;2.58349115827552,0.201788267359067,0.334671816614927;2.59349115827552,0.201010210851278,0.335967241279118;2.60349115827552,0.200238131367880,0.337262665943309;2.61349115827552,0.199471960299216,0.338558090607500;2.62349115827552,0.198711630081712,0.339853515271691;2.63349115827552,0.197957074178012,0.341148939935881;2.64349115827552,0.197208227057574,0.342444364600072;2.65349115827552,0.196465024177690,0.343739789264263;2.66349115827552,0.195727401964951,0.345035213928454;2.67349115827552,0.194995297797114,0.346330638592645;2.68349115827552,0.194268649985379,0.347626063256836;2.69349115827552,0.193547397757064,0.348921487921027;2.70349115827552,0.192831481238659,0.350216912585218;2.71349115827552,0.192120841439261,0.351512337249409;2.72349115827552,0.191415420234365,0.352807761913600;2.73349115827552,0.190715160350023,0.354103186577791;2.74349115827552,0.190020005347337,0.355398611241982;2.75349115827552,0.189329899607298,0.356694035906173;2.76349115827552,0.188644788315951,0.357989460570363;2.77349115827552,0.187964617449882,0.359284885234554;2.78349115827552,0.187289333762017,0.360580309898745;2.79349115827552,0.186618884767728,0.361875734562936;2.80349115827552,0.185953218731233,0.363171159227127;2.81349115827552,0.185292284652289,0.364466583891318;2.82349115827552,0.184636032253166,0.365762008555509;2.83349115827552,0.183984411965896,0.367057433219700;2.84349115827552,0.183337374919797,0.368352857883891;2.85349115827552,0.182694872929251,0.369648282548082;2.86349115827552,0.182056858481742,0.370943707212273;2.87349115827552,0.181423284726147,0.372239131876463;2.88349115827552,0.180794105461264,0.373534556540654;2.89349115827552,0.180169275124582,0.374829981204845;2.90349115827552,0.179548748781283,0.376125405869036;2.91349115827552,0.178932482113470,0.377420830533227;2.92349115827552,0.178320431409615,0.378716255197418;2.93349115827552,0.177712553554226,0.380011679861609;2.93700000000000,0.177439803992507,0.380336681406453;2.94700000000000,0.176237640759336,0.380336681406453;2.95700000000000,0.175047653339558,0.380336681406453;2.96700000000000,0.173869677860288,0.380336681406453;2.97700000000000,0.172703553196322,0.380336681406453;2.98700000000000,0.171549120915043,0.380336681406453;2.99700000000000,0.170406225222603,0.380336681406453;3.00700000000000,0.169274712911359,0.380336681406453;3.01700000000000,0.168154433308528,0.380336681406453;3.02700000000000,0.167045238226018,0.380336681406453;3.03700000000000,0.165946981911430,0.380336681406453;3.04700000000000,0.164859521000168,0.380336681406453;3.05700000000000,0.163782714468651,0.380336681406453;3.06700000000000,0.162716423588593,0.380336681406453;3.07700000000000,0.161660511882318,0.380336681406453;3.08700000000000,0.160614845079084,0.380336681406453;3.09700000000000,0.159579291072401,0.380336681406453;3.10700000000000,0.158553719878299,0.380336681406453;3.11700000000000,0.157538003594541,0.380336681406453;3.12700000000000,0.156532016360742,0.380336681406453;3.13700000000000,0.155535634319375,0.380336681406453;3.14700000000000,0.154548735577651,0.380336681406453;3.15700000000000,0.153571200170232,0.380336681406453;3.16700000000000,0.152602910022777,0.380336681406453;3.17700000000000,0.151643748916284,0.380336681406453;3.18700000000000,0.150693602452212,0.380336681406453;3.19700000000000,0.149752358018376,0.380336681406453;3.20700000000000,0.148819904755572,0.380336681406453;3.21700000000000,0.147896133524936,0.380336681406453;3.22700000000000,0.146980936876005,0.380336681406453;3.23700000000000,0.146074209015472,0.380336681406453;3.24700000000000,0.145175845776613,0.380336681406453;3.25700000000000,0.144285744589363,0.380336681406453;3.26700000000000,0.143403804451041,0.380336681406453;3.27700000000000,0.142529925897696,0.380336681406453;3.28700000000000,0.141664010976064,0.380336681406453;3.29700000000000,0.140805963216112,0.380336681406453;3.30700000000000,0.139955687604173,0.380336681406453;3.31700000000000,0.139113090556641,0.380336681406453;3.32700000000000,0.138278079894220,0.380336681406453;3.33700000000000,0.137450564816710,0.380336681406453;3.34700000000000,0.136630455878324,0.380336681406453;3.35700000000000,0.135817664963517,0.380336681406453;3.36700000000000,0.135012105263313,0.380336681406453;3.37700000000000,0.134213691252128,0.380336681406453;3.38700000000000,0.133422338665069,0.380336681406453;3.39700000000000,0.132637964475698,0.380336681406453;3.40700000000000,0.131860486874254,0.380336681406453;3.41700000000000,0.131089825246321,0.380336681406453;3.42700000000000,0.130325900151932,0.380336681406453;3.43700000000000,0.129568633305094,0.380336681406453;3.44700000000000,0.128817947553731,0.380336681406453;3.45700000000000,0.128073766860032,0.380336681406453;3.46700000000000,0.127336016281198,0.380336681406453;3.47700000000000,0.126604621950569,0.380336681406453;3.48700000000000,0.125879511059139,0.380336681406453;3.49700000000000,0.125160611837431,0.380336681406453;3.50700000000000,0.124447853537736,0.380336681406453;3.51700000000000,0.123741166416715,0.380336681406453;3.52700000000000,0.123040481718327,0.380336681406453;3.53700000000000,0.122345731657111,0.380336681406453;3.54700000000000,0.121656849401790,0.380336681406453;3.55700000000000,0.120973769059196,0.380336681406453;3.56700000000000,0.120296425658516,0.380336681406453;3.57700000000000,0.119624755135839,0.380336681406453;3.58700000000000,0.118958694319010,0.380336681406453;3.59700000000000,0.118298180912773,0.380336681406453;3.60700000000000,0.117643153484207,0.380336681406453;3.61700000000000,0.116993551448440,0.380336681406453;3.62700000000000,0.116349315054641,0.380336681406453;3.63700000000000,0.115710385372276,0.380336681406453;3.64700000000000,0.115076704277635,0.380336681406453;3.65700000000000,0.114448214440610,0.380336681406453;3.66700000000000,0.113824859311730,0.380336681406453;3.67700000000000,0.113206583109438,0.380336681406453;3.68700000000000,0.112593330807614,0.380336681406453;3.69700000000000,0.111985048123327,0.380336681406453;3.70700000000000,0.111381681504828,0.380336681406453;3.71700000000000,0.110783178119759,0.380336681406453;3.72700000000000,0.110189485843588,0.380336681406453;3.73700000000000,0.109600553248260,0.380336681406453;3.74700000000000,0.109016329591060,0.380336681406453;3.75700000000000,0.108436764803682,0.380336681406453;3.76700000000000,0.107861809481501,0.380336681406453;3.77700000000000,0.107291414873046,0.380336681406453;3.78700000000000,0.106725532869661,0.380336681406453;3.79700000000000,0.106164115995367,0.380336681406453;3.80700000000000,0.105607117396896,0.380336681406453;3.81700000000000,0.105054490833923,0.380336681406453;3.82700000000000,0.104506190669465,0.380336681406453;3.83700000000000,0.103962171860458,0.380336681406453;3.84700000000000,0.103422389948507,0.380336681406453;3.85700000000000,0.102886801050805,0.380336681406453;3.86700000000000,0.102355361851209,0.380336681406453;3.87700000000000,0.101828029591486,0.380336681406453;3.88700000000000,0.101304762062710,0.380336681406453;3.89700000000000,0.100785517596815,0.380336681406453;3.90700000000000,0.100270255058301,0.380336681406453;3.91700000000000,0.0997589338360881,0.380336681406453;3.92700000000000,0.0992515138355095,0.380336681406453;3.93700000000000,0.0987479554704557,0.380336681406453;3.94700000000000,0.0982482196556517,0.380336681406453;3.95700000000000,0.0977522677990728,0.380336681406453;3.96700000000000,0.0972600617944935,0.380336681406453;3.97700000000000,0.0967715640141680,0.380336681406453;3.98700000000000,0.0962867373016381,0.380336681406453;3.99700000000000,0.0958055449646682,0.380336681406453;4.00700000000000,0.0953279507683020,0.380336681406453;4.01700000000000,0.0948539189280413,0.380336681406453;4.02700000000000,0.0943834141031427,0.380336681406453;4.03700000000000,0.0939164013900308,0.380336681406453;4.04700000000000,0.0934528463158239,0.380336681406453;4.05700000000000,0.0929927148319735,0.380336681406453;4.06700000000000,0.0925359733080107,0.380336681406453;4.07700000000000,0.0920825885254019,0.380336681406453;4.08700000000000,0.0916325276715081,0.380336681406453;4.09700000000000,0.0911857583336483,0.380336681406453;4.10700000000000,0.0907422484932632,0.380336681406453;4.11700000000000,0.0903019665201789,0.380336681406453;4.12700000000000,0.0898648811669665,0.380336681406453;4.13700000000000,0.0894309615633985,0.380336681406453;4.14700000000000,0.0890001772109974,0.380336681406453;4.15700000000000,0.0885724979776769,0.380336681406453;4.16700000000000,0.0881478940924723,0.380336681406453;4.17700000000000,0.0877263361403596,0.380336681406453;4.18700000000000,0.0873077950571611,0.380336681406453;4.19700000000000,0.0868922421245353,0.380336681406453;4.20700000000000,0.0864796489650510,0.380336681406453;4.21700000000000,0.0860699875373421,0.380336681406453;4.22700000000000,0.0856632301313429,0.380336681406453;4.23700000000000,0.0852593493636030,0.380336681406453;4.24700000000000,0.0848583181726773,0.380336681406453;4.25700000000000,0.0844601098145935,0.380336681406453;4.26700000000000,0.0840646978583931,0.380336681406453;4.27700000000000,0.0836720561817454,0.380336681406453;4.28700000000000,0.0832821589666332,0.380336681406453;4.29700000000000,0.0828949806951083,0.380336681406453;4.30700000000000,0.0825104961451161,0.380336681406453;4.31700000000000,0.0821286803863884,0.380336681406453;4.32700000000000,0.0817495087764016,0.380336681406453;4.33700000000000,0.0813729569564005,0.380336681406453;4.34700000000000,0.0809990008474861,0.380336681406453;4.35700000000000,0.0806276166467666,0.380336681406453;4.36700000000000,0.0802587808235688,0.380336681406453;4.37700000000000,0.0798924701157114,0.380336681406453;4.38700000000000,0.0795286615258369,0.380336681406453;4.39700000000000,0.0791673323178018,0.380336681406453;4.40700000000000,0.0788084600131247,0.380336681406453;4.41700000000000,0.0784520223874900,0.380336681406453;4.42700000000000,0.0780979974673068,0.380336681406453;4.43700000000000,0.0777463635263226,0.380336681406453;4.44700000000000,0.0773970990822901,0.380336681406453;4.45700000000000,0.0770501828936854,0.380336681406453;4.46700000000000,0.0767055939564790,0.380336681406453;4.47700000000000,0.0763633115009559,0.380336681406453;4.48700000000000,0.0760233149885864,0.380336681406453;4.49700000000000,0.0756855841089446,0.380336681406453;4.50700000000000,0.0753500987766753,0.380336681406453;4.51700000000000,0.0750168391285074,0.380336681406453;4.52700000000000,0.0746857855203142,0.380336681406453;4.53700000000000,0.0743569185242183,0.380336681406453;4.54700000000000,0.0740302189257407,0.380336681406453;4.55700000000000,0.0737056677209950,0.380336681406453;4.56700000000000,0.0733832461139231,0.380336681406453;4.57700000000000,0.0730629355135742,0.380336681406453;4.58700000000000,0.0727447175314244,0.380336681406453;4.59700000000000,0.0724285739787378,0.380336681406453;4.60700000000000,0.0721144868639673,0.380336681406453;4.61700000000000,0.0718024383901948,0.380336681406453;4.62700000000000,0.0714924109526102,0.380336681406453;4.63700000000000,0.0711843871360279,0.380336681406453;4.64700000000000,0.0708783497124418,0.380336681406453;4.65700000000000,0.0705742816386158,0.380336681406453;4.66700000000000,0.0702721660537109,0.380336681406453;4.67700000000000,0.0699719862769479,0.380336681406453;4.68700000000000,0.0696737258053052,0.380336681406453;4.69700000000000,0.0693773683112497,0.380336681406453;4.70700000000000,0.0690828976405032,0.380336681406453;4.71700000000000,0.0687902978098400,0.380336681406453;4.72700000000000,0.0684995530049187,0.380336681406453;4.73700000000000,0.0682106475781453,0.380336681406453;4.74700000000000,0.0679235660465674,0.380336681406453;4.75700000000000,0.0676382930897999,0.380336681406453;4.76700000000000,0.0673548135479814,0.380336681406453;4.77700000000000,0.0670731124197594,0.380336681406453;4.78700000000000,0.0667931748603064,0.380336681406453;4.79700000000000,0.0665149861793635,0.380336681406453;4.80700000000000,0.0662385318393135,0.380336681406453;4.81700000000000,0.0659637974532814,0.380336681406453;4.82700000000000,0.0656907687832630,0.380336681406453;4.83700000000000,0.0654194317382794,0.380336681406453;4.84700000000000,0.0651497723725597,0.380336681406453;4.85700000000000,0.0648817768837483,0.380336681406453;4.86700000000000,0.0646154316111396,0.380336681406453;4.87700000000000,0.0643507230339362,0.380336681406453;4.88700000000000,0.0640876377695335,0.380336681406453;4.89700000000000,0.0638261625718281,0.380336681406453;4.90700000000000,0.0635662843295504,0.380336681406453;4.91700000000000,0.0633079900646210,0.380336681406453;4.92700000000000,0.0630512669305304,0.380336681406453;4.93700000000000,0.0627961022107415,0.380336681406453;4.94700000000000,0.0625424833171151,0.380336681406453;4.95700000000000,0.0622903977883571,0.380336681406453;4.96700000000000,0.0620398332884878,0.380336681406453;4.97700000000000,0.0617907776053327,0.380336681406453;4.98700000000000,0.0615432186490347,0.380336681406453;4.99700000000000,0.0612971444505865,0.380336681406453]}};

    if options.Controls.DET==1
        switch self.isolationType
        case 'LRB'

            microFieldsPar{4} = {'costRatio','DS','DLR','Mu','HeightBounds','ARatioBounds','SigmaBounds','NseedsHeight','NseedsARatio','NseedsSigma','Diameter_height','Glead','Grubber','epsYlead','fylead'};
            microFieldsParVals{4} = {0.02,[0.4 1 1.25],[0 20 100 100]/100,13.33,[0.3,0.8],[10,30],[4,10],10,10,10,2,130,1,0.075,10};

        case 'FPS'

            microFieldsPar{4} = {'costRatio','DS','DLR','Dyiso','MiuBounds','RadiusBounds','DiameterBounds','NseedsMiu','NseedsDiameter','NseedsRadius'};
            microFieldsParVals{4} = {0.02,[0.4 1 1.25],[0 20 100 100]/100,0.001,[0.02,0.15],[3,20],[0.6,1.6],10,10,10};

        otherwise,  warning('IsolationType not recognized, input LRB or FPS');  

        end
    end

    if options.Controls.DET==0

        switch self.isolationType
            
            case 'LRB'
            microFieldsPar{4} = {'costRatio','DS','DLR','Mu','fyBounds','k1Bounds','alphaBounds','NseedsFy','NseedsK1','NseedsAlpha','Diameter_height','Glead','Grubber','epsYlead','fylead'};
            microFieldsParVals{4} = {0.002,[0.4 1 1.25],13.33,[0 20 100 100]/100,[0.05,0.15],[0.5,15],[0.03,0.22],10,10,10,2,130,1,0.075,10};

            case 'FPS'
            microFieldsPar{4} = {'costRatio','DS','DLR','Dyiso','fyBounds','alphaBounds','DiameterBounds','NseedsFy','NseedsAlpha','NseedsDiameter'};
            microFieldsParVals{4} = {0.002,[0.4 1 1.25],[0 20 100 100]/100,0.002,[0.005,0.12],[0.01,0.3],[0.3,0.8],10,10,10};

            otherwise,  warning('IsolationType not recognized, input LRB or FPS') 

        end
    end
    
    switch self.ssType
        case 'Frame'
        microFieldsPar{5} = {'costRatio','DS','DLR','FyBounds','NseedsFy','numStoreys','FloorsWeight','Gamma','InterStoreyHeights','Material','epsYsteel','lengthBeams','heightBeams','t1','DYatTop','DYeff','partFactor','meffFactor','drifty','accProf'};
        microFieldsParVals{5} = {0.124,[0.33 0.66 1],[0 5 10 15]/100,[0.1,0.4],3,10,[80,100,100,100],0.75,[4,3,3],'Steel',0.002,6,0.3,0.4,0,0,1,1,0,'Modal'};


        case 'Wall'  
        microFieldsPar{5} = {'costRatio','DS','DLR','FyBounds','NseedsFy','numStoreys','FloorsWeight','Gamma','InterStoreyHeights','epsYsteel','lengthWall','t1','DYatTop','DYeff','partFactor','meffFactor','drifty','accProf'};
        microFieldsParVals{5} = {0.124,[0.33 0.66 1],[0 5 10 15]/100,[0.1,0.4],3,10,[80,100,100,100],0.75,[4,3,3],0.002,4,0.3,0,0,01,1,0,'Modal'};

        otherwise, warning('Variable ssType not recognized, input frame or wall');
    end
    
    microFieldsPar{6} = {'explicitNSC', 'Dsens', 'Asens'};
    microFieldsParVals{6} = {false, struct, struct};
    
    microFieldsPar{7} = {'DET','CalcHeffstar','ColLosses'};
    microFieldsParVals{7} = {1,1,0};

    microFieldsPar{8} = {'replacementCost','indLossesPerDay'};
    microFieldsParVals{8} = {1,100000};

    microFieldsPar{9} = {'numWorkersPerFloor','workingHoursRatio','workDaysRatio','repairScheme'};
    microFieldsParVals{9} = {12,0.33,0.71,'series'};

    for F = 1 : numel(macroFieldsPar)
        for f = 1 : numel(microFieldsPar{F})
            self.parameters.(macroFieldsPar{F}).(microFieldsPar{F}{f}) = ...
                microFieldsParVals{F}{f};
        end
    end

    % overwrite fields if some parameter is specified
    macroFieldsOptional = fieldnames(options);
    for OF = 1 : numel(macroFieldsOptional)
        microFieldsOptional = fieldnames(options.(macroFieldsOptional{OF}));
        for of = 1 : numel(microFieldsOptional)
            if isfield(self.parameters.(macroFieldsOptional{OF}), microFieldsOptional{of}) == 1
                self.parameters.(macroFieldsOptional{OF}).(microFieldsOptional{of}) = options.(macroFieldsOptional{OF}).(microFieldsOptional{of});
            else
                error('Field %s.%s is not allowed', macroFieldsOptional{OF}, microFieldsOptional{of})
            end
        end
    end
end

function self = setSeedSDoFproperties(self)

    Iso=self.parameters.Iso;

    self.fySS=linspace(self.parameters.SS.FyBounds(1),self.parameters.SS.FyBounds(2),self.parameters.SS.NseedsFy)';


   if strcmp(self.isolationType,'LRB')

        self.parameters.Iso=Iso;

        fyDummy=linspace(self.parameters.Iso.fyBounds(1),self.parameters.Iso.fyBounds(2),self.parameters.Iso.NseedsFy);
        k1Dummy=linspace(self.parameters.Iso.k1Bounds(1),self.parameters.Iso.k1Bounds(2),self.parameters.Iso.NseedsK1);
        alphaDummy=linspace(self.parameters.Iso.alphaBounds(1),self.parameters.Iso.alphaBounds(2),self.parameters.Iso.NseedsAlpha);

        [FYISO, K1ISO,ALPHAISO]= meshgrid(fyDummy,k1Dummy,alphaDummy);


        self.fyiso = reshape(FYISO, numel(FYISO), 1);
        self.k1 = reshape(K1ISO, numel(K1ISO), 1);
        self.alphaiso = reshape(ALPHAISO, numel(ALPHAISO), 1);


        self.t1=2*pi*(1./(9.81*self.k1)).^0.5;
        self.t2 = 2*pi*((1/9.81.*ones(length(self.fyiso),1))./(self.k1.*self.alphaiso)).^0.5;
        self.k2= self.k1.*self.alphaiso;
        % Simplified calculation of parameters from this line on to
        % be corrected in detailing of LRB
        self.Dyiso = self.fyiso./self.k1; % to be corrected if flat sliders are used (fy=fy+miu*Wsliders/Wtot)
        self.hiso=Iso.Glead*(self.fyiso/Iso.fylead)./(self.k1-self.k1.*self.alphaiso);
        self.muiso=self.hiso./self.Dyiso;
        self.DuIso=self.Dyiso.*self.muiso; % Simplified assumption
        
        %

    elseif strcmp(self.isolationType,'FPS') 

        fyDummy=linspace(self.parameters.Iso.fyBounds(1),self.parameters.Iso.fyBounds(2),self.parameters.Iso.NseedsFy);
        alphaDummy=linspace(self.parameters.Iso.alphaBounds(1),self.parameters.Iso.alphaBounds(2),self.parameters.Iso.NseedsAlpha);
        DuIsoDummy=linspace(self.parameters.Iso.deltaUBounds(1),self.parameters.Iso.deltaUBounds(2),self.parameters.Iso.NseedsDeltaU);

        [self.FYISO, self.ALPHAISO, self.DU]= meshgrid(fyDummy,alphaDummy,DuIsoDummy);

        self.fyiso = reshape(self.FYISO, numel(self.FYISO), 1);
        self.alphaiso = reshape(self.ALPHAISO, numel(self.ALPHAISO), 1);            
        self.DuIso = reshape(self.DU, numel(self.DU), 1);

        self.Dyiso =cell(numel(self.fyiso),1); 
        self.Dyiso(:,1) = {Iso.Dyiso}; 

        self.k1=self.fyiso./Iso.Dyiso;
        self.t1=2*pi*((1/9.81.*ones(length(self.fyiso),1))./(self.k1)).^0.5;
        self.t2 = 2*pi*((1/9.81.*ones(length(self.fyiso),1))./(self.k1.*self.alphaiso)).^0.5;
        self.k2= self.k1.*self.alphaiso;                
        self.muiso=self.DuIso./self.Dyiso;

    end


        self.hyst = cell(numel(self.fyiso)); 
        self.hyst(:) = {self.parameters.General.hysteresis};
        self = checkExtrapolations(self);

end

function self = setSeedSDoFproperties2(self)

    Iso=self.parameters.Iso;
    SS=self.parameters.SS;

    self.fySS=linspace(SS.FyBounds(1),SS.FyBounds(2),SS.NseedsFy)';


    if strcmp(self.isolationType,'LRB')

        heightDummy=linspace(Iso.HeightBounds(1),Iso.HeightBounds(2),Iso.NseedsHeight);
        aratioDummy=linspace(Iso.ARatioBounds(1),Iso.ARatioBounds(2),Iso.NseedsARatio);
        sigmaDummy= linspace(Iso.SigmaBounds(1),Iso.SigmaBounds(2),Iso.NseedsSigma);

        COMB=combvec(heightDummy,aratioDummy,sigmaDummy);

        self.hiso=COMB(1,:)';
        self.aratio=COMB(2,:)';
        self.sigma=COMB(3,:)';

        % Refine ultimate displacement calculation for LRB in LRBgetparameters

        Wtot=sum(self.parameters.SS.FloorsWeight);

        for i=1:length(COMB)
            [self.Dyiso(i,1),self.k1(i,1),self.alphaiso(i,1),self.fyiso(i,1),self.muiso(i,1),self.areaiso(i,1),self.numiso(i,1),self.Teff(i,1),self.Deff(i,1)]=...
                self.LRBgetparameters(Iso.Diameter_height,self.hiso(i),Iso.Grubber,Iso.Glead,self.aratio(i),Iso.epsYlead,self.sigma(i),Wtot);
        end
        
        self.DuIso=self.muiso.*self.Dyiso;

%                 [self.FYISO, self.K1ISO, self.ALPHAISO]= meshgrid(self.fyiso,self.alphaiso,self.alphaiso);

    elseif strcmp(self.isolationType,'FPS') 

        fyDummy=linspace(Iso.MiuBounds(1),Iso.MiuBounds(2),Iso.NseedsMiu);
        risoDummy=linspace(Iso.RadiusBounds(1),Iso.RadiusBounds(2),Iso.NseedsRadius);
        diamDummy=linspace(Iso.DiameterBounds(1),Iso.DiameterBounds(2),Iso.NseedsDiameter);

        COMB=combvec(fyDummy,risoDummy,diamDummy);

        self.fyiso=COMB(1,:)';
        self.Riso=COMB(2,:)';
        self.diamFPS=COMB(3,:)';

        for i=length(COMB):-1:1
            [self.k1(i,1),self.alphaiso(i,1),self.muiso(i,1)]=...
                self.FPSgetparameters(self.fyiso(i),self.Riso(i),self.diamFPS(i),Iso.Dyiso);
        end            
    self.Dyiso=Iso.Dyiso*ones(length(self.fyiso),1);
    self.DuIso=self.diamFPS./2;
    end
    

    self.t1=2*pi*((1/9.81.*ones(length(self.fyiso),1))./(self.k1)).^0.5;
    self.t2 = 2*pi*((1/9.81.*ones(length(self.fyiso),1))./(self.k1.*self.alphaiso)).^0.5;
    self.k2= self.k1.*self.alphaiso;
    self.hyst = cell(numel(self.fyiso)); 
    self.hyst(:) = {self.parameters.General.hysteresis};
    self = checkExtrapolations(self);
    self.DncIso=self.DuIso.*self.parameters.Iso.DS(end);

end

function self = calcIsolationEff(self)

    SS=self.parameters.SS;
   






end

function self = setIsoDSThreshold(self)

    for ds = numel(self.parameters.Iso.DS) : -1 : 1

        self.ductThresholdsIso(:,ds) = ...
            self.parameters.Iso.DS(ds) .* self.muiso;

    end 
end


function self = calcHeffstar(self)
%Calculate the effective height of the isolated structure for
%each IM level considering the median peak displacement and...
%acceleration of the isolation layer. 
    if self.parameters.Controls.CalcHeffstar==1

    SS=self.parameters.SS;
    tol=0.05;
    diff=1;


    % Get Isolation layer displacement 
    for n = numel(self.fyiso):-1:1
        EC = self.IMdef./self.fyiso(n) <= 1; %ElasticCases
        dispIso(EC,1)=(self.IMdef(EC)./self.fyiso(n))*...
            self.Dyiso(n);
        dispIso(~EC,1)=(1+self.powerLawDuct(n,1).*...
            ((self.IMdef(~EC)./self.fyiso(n))-1)).*...
            self.Dyiso(n);
        %Control plot
        %figure()
        %plot(self.IMdef./self.fyiso(n),dispIso)
        accIso(EC,1)=self.IMdef(EC);
        accIso(~EC,1)=(1+self.powerLawAccRat(n,1).*...
            ((self.IMdef(~EC)./self.fyiso(n))-1)).*...
            self.fyiso(n); 
        % Control plot
        %figure()
        %plot(self.IMdef,accIso)
        
        for SSt=numel(self.fySS):-1:1
            heffZero=SS.heff;
            i=2;
            for m= numel(self.IMdef):-1:1    
                heffisolated=heffZero+tol;
                i=2;
                while abs(diff)>tol && i<5
                     kstar=self.fySS(SSt)*(SS.heff+SS.StoreyHeights(1))/...
                       (SS.dY_func(heffisolated(i-1))*heffisolated(i-1));
%                     kstar=self.fySS(SSt)*SS.Gamma/SS.dY_func(heffisolated(i-1));
                    DHeffstar=SS.Gamma*accIso(m,1)/(kstar);
                    fractionDy=DHeffstar/(SS.dY_func(heffisolated(i-1)));
                    SSdisp=fractionDy.*SS.dY_func(SS.StoreyHeights(2:end));
                    heffisolated(i)=...
                        sum(SS.FloorsWeight.*SS.StoreyHeights.*[dispIso(m),SSdisp+dispIso(m)])/...
                        sum(SS.FloorsWeight.*[dispIso(m),SSdisp+dispIso(m)]);
                    diff=heffisolated(i-1)-heffisolated(i);
                    i=i+1;
                    conv(m,SSt)= abs(diff)<tol;
                end 
                self.expDisp(m,:)=[dispIso(m),SSdisp+dispIso(m)];
                self.expVss(m,n,SSt)= fractionDy*self.fySS(SSt);
                diff=1;
                self.heffstar(m,n,SSt)=heffisolated(i-1)-SS.StoreyHeights(1);
                heffZero=heffisolated(i-1);
                clear heffisolated
            end
        end
    end

    else
        self.heffstar=ones(numel(self.IMdef),numel(self.fyiso),...
            numel(self.fySS)).*self.parameters.SS.heff;
    end
end

function self=ductSSdonatello(self)
    SS=self.parameters.SS;
    for SSt=numel(self.fySS):-1:1
        for n = numel(self.fyiso):-1:1 
            R=self.IMdef./self.fyiso(n);
            EC = R <= 1; %ElasticCases
            self.PSDMSS(EC,n,SSt)=(R(EC).*self.fyiso(n)*9.81*SS.t1(1)^2*SS.meffFactor/...
                (4*pi^2*SS.partFactor))./...
                SS.DYeff;
            self.PSDMSS(~EC,n,SSt)=(self.fyiso(n).*...
                (1+self.powerLawAccRat(n,1)*(R(~EC)-1))*9.81*SS.t1(1)^2*SS.meffFactor/...
                (4*pi^2*SS.partFactor))./...
                SS.DYeff;
        end
    end

end

function self = ductSSYe(self)

        SS=self.parameters.SS;

    for SSt=numel(self.fySS):-1:1
        for n = numel(self.fyiso):-1:1 
            R=self.IMdef./self.fyiso(n);
            EC = R <= 1; %ElasticCases
            Q=self.fyiso(n)-(self.k2(n)*self.Dyiso(n));
            Diso(EC)=R(EC).*self.Dyiso(n);
            Diso(~EC)=self.Dyiso(n).*(1+self.powerLawDuct(n,1)*(R(~EC)-1));
            Teff(EC)=self.t1(n);
            keff(~EC)=self.k2(n)+Q./Diso(~EC);
            Teff(~EC)=2*pi.*(1./(keff(~EC)*9.81)).^0.5;
            Trat=Teff./SS.t1(SSt);
            self.PSDMSS(:,n,SSt)=Diso./...
                (Trat.^2).*sum(SS.FloorsWeight)/...
                (SS.FloorsWeight(1)+sum(SS.FloorsWeight(2:end))*SS.meffFactor)./...
                SS.DYeff;
        end
    end
end



function self=calcSSPSDM(self)
    SS=self.parameters.SS;
    for SSt=numel(self.fySS):-1:1
        for n = numel(self.fyiso):-1:1    
            R=self.IMdef./self.fyiso(n);
            EC = R <= 1; %ElasticCases
%             Vpf = SS.FloorsWeight.*self.expDisp(:,:)./...
%                 (sum(SS.FloorsWeight.*self.expDisp(:,:)));
%             Vss_Vtot = sum(Vpf(:,2:end),2)./sum(Vpf,2);
            self.PSDMSS(EC,n,SSt)=R(EC).*self.fyiso(n)./...
                (self.fySS(SSt)*SS.heff./self.heffstar(EC,n,SSt)*SS.partFactor);
            self.PSDMSS(~EC,n,SSt)=self.fyiso(n).*...
                (1+self.powerLawAccRat(n,1)*(R(~EC)-1))./...
                (self.fySS(SSt).*SS.heff./self.heffstar(~EC,n,SSt)*SS.partFactor);
        end
    end
end

function self = checkExtrapolations(self)

    maxPeriodSeeds = max(max(self.t1));
    minPeriodSeeds = min(min(self.t1));

    maxPeriodHazModel = max(self.parameters.Hazard.periodsHazard);
    minPeriodHazModel = min(self.parameters.Hazard.periodsHazard);

    if maxPeriodSeeds > maxPeriodHazModel
        warning(['Some hazard curve is extrapolated\n' ...
            'Max isolation period Seeds: %2.2fs \n' ...
            'Max period Hazard model: %2.2fs'], ...
            maxPeriodSeeds, maxPeriodHazModel)
    end

    if minPeriodSeeds < minPeriodHazModel
        warning(['Some hazard curve is extrapolated\n' ...
            'Min isolation period Seeds: %2.2fs \n' ...
            'Min period Hazard model: %2.2fs'], ...
            minPeriodSeeds, minPeriodHazModel)
    end

    % To be added check for extrapolation of Fy, T1 and alpha outside...
    % the training dataset boundaries
end   

function self = plotCandidatesScatter(self,C,title)
       
    figure('Position', [10,10,1000,500]);
    hold on
    if length(self.fySS)==1;numCol=1;else; numCol=length(self.fySS)/2; end
    tiledlayout(numCol,2)
    

        for SS=1:length(self.fySS)
            switch self.isolationType
            case 'LRB'
            xlabel='height_{iso} (m)';
            ylabel='F_y (g)';
            zlabel='h_{iso}';
            clabel='EAL';
    
            x=self.hiso(C(:,SS)==1);
            y=self.fyiso(C(:,SS)==1);
            z=self.alphaiso(C(:,SS)==1);
            c=self.eal(C(:,SS)==1,SS);
            xbounds=[min(self.hiso) max(self.hiso)];
            ybounds=[0.0 max(self.fyiso)];
            zbounds=[0.0 max(self.alphaiso)];
    
            case 'FPS'
    
            xlabel='Diameter (m)';
            ylabel='miu (g)';
            zlabel='Radius (m)';
            clabel='EAL NSC';
            x=self.diamFPS(C(:,SS)==1);
            y=self.fyiso(C(:,SS)==1);
            z=self.Riso(C(:,SS)==1);
            c=self.eal(C(:,SS)==1,SS);
            xbounds=[min(self.diamFPS) max(self.diamFPS)];
            ybounds=[0.0 max(self.fyiso)];
            zbounds=[0.0 max(self.Riso)];
    
            end

            Title=strcat(title,' for fy_{SS}=',...
                num2str(round(self.fySS(SS),2)),' g');
            nexttile;
            grid on
            lossBasedDesignISO.plot4Dscattered...
                (x,y,z,c,xlabel,ylabel,zlabel,clabel ,Title,...
                xbounds,ybounds,zbounds,...
                'auto','normal','reverse','normal');

        hold off
        end
end


function self = candidatesInfo(self)
   
    switch self.isolationType

        case 'LRB'
            
            headers = {'indIso','indSS','EAL','Hiso','Arubber','Alead','Aratio','Sigma'};

            for SS=length(self.fySS):-1:1
                indCAND=self.indCANDIDATES{SS};
                indSS = ones(length(indCAND),1).*SS;
                HisoCand=self.hiso(indCAND); 
                Arubber=self.det.Arubber(indCAND);
                Alead=self.det.Alead(indCAND);
                ARatioCand=self.aratio(indCAND);
                SigmaCand=self.sigma(indCAND);
                EalCand=self.eal(indCAND);
                self.designpar{SS}=[indCAND, indSS, EalCand,HisoCand, Arubber, Alead ARatioCand SigmaCand];
                T = array2table(self.designpar{SS}, 'VariableNames', headers);
                disp(T)
            end

            

        case 'FPS'
        
            headers = {'indIso','indSS','Diameter','alpha','Fyiso'};

            for SS=length(self.fySS):-1:1
                indCAND=self.indCANDIDATES{SS};
                indSS = ones(length(indCAND),1).*SS;
                DiamCand=self.diamFPS(indCAND); 
                Fyiso=self.fyiso(indCAND);
                alpha=self.alphaiso(indCAND);
                self.designpar{SS}=[indCAND, indSS ...
                     DiamCand alpha Fyiso];
                T = array2table(self.designpar{SS}, 'VariableNames', headers);
                disp(T)
            end
    end

    

end

function printDesigninfoLRB(self)
    indIso=self.Selected(1);
    indSS=self.Selected(2);
    % Properties table
    headers={'indIso','indSS','hiso','#iso','t1','fyiso','k1','alphaiso',...
        'fySS', 'DyTopSS'};
    row=[indIso indSS ...
        self.hiso(indIso) self.numiso(indIso) self.t1(indIso) ...
        self.fyiso(indIso) self.k1(indIso) self.alphaiso(indIso) self.fySS(indSS)...
        self.parameters.SS.Dy];
    SelectionInfo=array2table(row, 'VariableNames', headers)
    % Results table
   headers={'indIso','indSS','EAL','EALiso','EALSS','EALNSCD','EALNSCA','EALCollapse',...
        'MAFEyieldSS','MAFEcol','isCAND'};
    row=[indIso indSS ...
        self.eal(indIso,indSS) self.ealIso(indIso) ...
        self.ealSS(indIso,indSS) self.ealNSCD(indIso,indSS)...
        self.ealNSCA(indIso,indSS) self.ealCollapse(indIso)...
        self.mafeDSyieldSS(indIso,indSS) ...
        self.mafeDSColIso(indIso) self.isCANDIDATE(indIso,indSS)]; 
    SelectionInfo=array2table(row, 'VariableNames', headers)
end



function printDesigninfoFPS(self)
    indIso=self.Selected(1);
    indSS=self.Selected(2);
    % Properties table
    headers={'indIso','indSS','Diameter','t1','fyiso','alphaiso',...
        'fySS', 'DyTopSS'};
    row=[indIso indSS ...
        self.diamFPS(indIso) self.t1(indIso) ...
        self.fyiso(indIso) self.alphaiso(indIso) self.fySS(indSS)...
        self.parameters.SS.Dy];
    SelectionInfo=array2table(row, 'VariableNames', headers)
    % Results table
    headers={'indIso','indSS','EAL','EALiso','EALSS','EALNSCD','EALNSCA',...
        'MAFEyieldSS','MAFEcol','isCAND'};
    row=[indIso indSS ...
        self.eal(indIso,indSS) self.ealIso(indIso) ...
        self.ealSS(indIso,indSS) self.ealNSCD(indIso,indSS)...
        self.ealNSCA(indIso,indSS) self.mafeDSyieldSS(indIso,indSS) ...
        self.mafeDSColIso(indIso) self.isCANDIDATE(indIso,indSS)]; 
    SelectionInfo=array2table(row, 'VariableNames', headers)
end

end

methods(Static)

function [Dy,k1_W,k2_k1,Fy_W,mu,A,numiso,Teff,Deff] = LRBgetparameters(...
        D_h,h,G_R,G_L,ARatios,ey_L,sigma,Wtot)

    D=h*D_h; %Design decision

    A=Wtot/(sigma*1000);
    numiso=A/(pi*D^2/4);
    Ap=A/(ARatios+1);
    Ar=A/(1+1/ARatios);

    Kr=G_R*1000*Ar/h;
    kp=G_L*1000*Ap/h;
    k1_W=(Kr+kp)/Wtot; 
    k2_W=Kr/Wtot;
    k2_k1=k2_W/k1_W;

    Dy=ey_L*h;
    Fy_W=k1_W*Dy;

    Du=h; % Need to refine the calculation
    mu=Du/Dy;

    Q=(Fy_W-Dy*k2_W)*Wtot;
    keff=Kr+Q/Du;
    Teff=2*pi*(Wtot/(9.81*keff))^0.5;
    Deff=4*Q*(Du-Dy)/(2*pi*(k2_W*Du*Wtot+Q)*Du);

end


function [k1_W,k2_k1,mu] = FPSgetparameters(...
        miu,R,DuIso,Dyiso)

     Fy_W=miu;

     k2_W=1/R;
     k1_W=Fy_W/Dyiso;
     k2_k1=k2_W/k1_W;

     mu=(DuIso/2)/Dyiso; %Approximation 

%              Q=Fy_W-k2_W*Dyiso;
%              T=2*pi*(R/9810)^0.5;            
%              Deq=2*miu/(pi*D/R+miu);


end      

function hazardExtrap = logExtrapHazardCurve(...
        hazardCurve, faultRate, mode)

    if nargin < 3; mode = 'linear'; end

    if strcmp(mode, 'powerLaw')
        % out of sensitivity these seem reasonable
        IMzero = 0.01;
        IMstep = IMzero;

        % define a (log)line between IM=0 and IM=IMmin
        x1 = log(IMzero);
        x2 = log(hazardCurve(1,1));
        y1 = log(faultRate);
        y2 = log(hazardCurve(1,2));

        k1 = (y2-y1)/(x2-x1);
        k2 = exp( -(y2-y1)/(x2-x1) * x1 + y1 );

        % add the points to the hazard curve
        imExtrap = IMzero : IMstep : hazardCurve(1,1);
        mafeExtrap = k2*imExtrap.^k1;

        hazardExtrap(:,1) = ...
            [0; imExtrap'; hazardCurve(:,1)];
        hazardExtrap(:,2) = ...
            [mafeExtrap(1); mafeExtrap'; hazardCurve(:,2)];
    else
        hazardExtrap(:,1) = [0; hazardCurve(:,1)];
        hazardExtrap(:,2) = [faultRate; hazardCurve(:,2)];
    end
end


function MAFEds = calculateMAFEds(...
        fragMedians, fragStdev, hazardCurve)

    IMdef = linspace(hazardCurve(1,1), hazardCurve(end,1), 1000);
    hazardCurveResampled = interp1(...
        hazardCurve(:,1), hazardCurve(:,2), IMdef)';

    for ds = numel(fragMedians) : -1 : 1
        fragilities(:,ds) = logncdf(IMdef,...
            log(fragMedians(ds)), fragStdev(ds));

        MAFEds(ds) = trapz(IMdef, ...
            fragilities(:,ds) .* hazardCurveResampled);
    end

end


function plot4Dscattered(x,y,z,c,xname,yname,zname,cname,Title,Xlim,Ylim,Zlim,Clim,xdir,ydir,zdir)

    scatter3(x,y,z,42,c,'filled')    
    ax = gca;
    ax.XDir = xdir;
    ax.YDir = ydir;
    ax.ZDir = zdir;
    view(-31,14)  
    xlabel(xname)
    ylabel(yname)
    zlabel(zname)
%     cb = colorbar;  
    cb = colorbar('southoutside');                                    
    cb.Label.String = cname;
    xlim(Xlim);
    ylim(Ylim);
    zlim(Zlim);
    clim(Clim);
    title(Title);
    set(gca, 'FontSize', 16);

end

function self = plotCandidatesBackbone(self)

    C= self.isISOCOLACC & self.isSSYIELDACC;   
    ind = self.Selected;
    colors=cool(3);
    %[1 0 0 ;
     %       0 0 1 ,
      %      0 1 0 ];

    figure(); hold on

    for n=length(self.fyiso):-1:1
            plot([0,self.Dyiso(n),self.DuIso(n)],[0,self.fyiso(n),...
            self.fyiso(n)+self.k2(n).*(self.DuIso(n)-self.Dyiso(n))],...
            'color',[0.7 0.7 0.7],'LineWidth',0.8)

    end

%     
     for SS=length(self.fySS):-1:1
        N=sum(C(:,SS));
%         nexttile;hold on
        grid on
        xlabel('Displacement [m]')
        ylabel('Fy_{iso} [-]')
        set(gca, 'FontSize', 12) 
        
        for n=length(self.fyiso):-1:1
            plot([0,self.Dyiso(n),self.DuIso(n)].*C(n,SS),[0,self.fyiso(n),...
            self.fyiso(n)+self.k2(n).*(self.DuIso(n)-self.Dyiso(n))].*C(n,SS),...
            'color',colors(SS,:),'LineWidth',1)
        end
     end

%      for n=length(self.fyiso):-1:1
%             scatter(self.DuIso(n),...
%             (self.fyiso(n)+self.k2(n).*(self.DuIso(n)-self.Dyiso(n))),...
%             80,[0.8 0.8 0.8],'filled','x','MarkerEdgeColor',[0.8 0.8 0.8]); 
%     end

      for SS=length(self.fySS):-1:1
          for n=length(self.fyiso):-1:1
              scatter(self.DuIso(n).*C(n,SS),...
        (self.fyiso(n)+self.k2(n).*(self.DuIso(n)-self.Dyiso(n))).*C(n,SS),...
         80,colors(SS,:),'filled','x','MarkerEdgeColor',[0 0 0]); 
          end
      end

%      for SS=[2]%length(self.fySS):-1:1
%          for ind=1:length(self.fyiso)
%     
%             plot([0,self.Dyiso(ind),self.DuIso(ind)].*~C(ind,SS),[0,self.fyiso(ind),...
%                 self.fyiso(ind)+self.k2(ind).*(self.DuIso(ind)-self.Dyiso(ind))].*~C(ind,SS),...
%                 '-','color',[0.6 0.6 0.6])
%          end
%      end




       legend('Reliability not complying',...
           'Reliability complying',...
           'Max. disp. capacity','Location','southeast');
end


function self = plotLossSurf(self,nc,yss,xData,yData,zData,ctData)
       
%   figure('Position', [10,10,1000,500]);
figure()
    hold on
%     if length(self.fySS)==1;numCol=1;else; numCol=length(self.fySS)/2; end
%     tiledlayout(numCol,2)
%  

        for SS=1%:length(self.fySS)
        
%         nexttile
%         Title=strcat('Loss surfaces for different -',{' '}, ctData,'.',{' '}, 'Fy_{SS} = ',...
%             num2str(round(self.fySS(SS),2)),' g');
%         title(Title)
        

        view(3)

        colors=cool((numel(unique(self.(ctData))))/1);
        i=1;

        uniqueCt=unique(self.(ctData));

        FaceAlpha=[0.6 0.6 1 0.6 0.6 0.6 0.6 0.6];

        for D=[1 2 3 4 8] %D=1:1:(numel(unique(self.(ctData))))

            d=uniqueCt(D);
            ind = self.(ctData)==d;

            NC=nc(ind);
            ySS=yss(ind,SS);

            x1=self.(xData)(ind);
            y1=self.(yData)(ind);
            z1=self.(zData)(ind,SS);

            x2=x1;
            y2=y1;
            z2=z1;
       
            x3=x1;
            y3=y1;
            z3=z1;            

            x4=x1;
            y4=y1;
            z4=z1;

            x2(~NC)=NaN;
            y2(~NC)=NaN;
            z2(~NC)=NaN;

            x3(~ySS(:))=NaN;
            y3(~ySS(:))=NaN;
            z3(~ySS(:))=NaN;

            x4(~NC|~ySS(:))=NaN;
            y4(~NC|~ySS(:))=NaN;
            z4(~NC|~ySS(:))=NaN;

            X1=reshape(x1, numel(unique(self.(xData))), []);
            Y1=reshape(y1, numel(unique(self.(yData))),[]);
            Z1=reshape(z1, numel(unique(self.(yData))), []);

            X2=reshape(x2, numel(unique(self.(xData))), []);
            Y2=reshape(y2, numel(unique(self.(yData))),[]);
            Z2=reshape(z2, numel(unique(self.(yData))), []);

            X3=reshape(x3, numel(unique(self.(xData))), []);
            Y3=reshape(y3, numel(unique(self.(yData))),[]);
            Z3=reshape(z3, numel(unique(self.(yData))), []);

            X4=reshape(x4, numel(unique(self.(xData))), []);
            Y4=reshape(y4, numel(unique(self.(yData))),[]);
            Z4=reshape(z4, numel(unique(self.(yData))), []);

            
            S1=surface(X1(1:2:end,1:2:end),Y1(1:2:end,1:2:end),Z1(1:2:end,1:2:end),'FaceAlpha',0,'EdgeAlpha',0,'EdgeColor',[0,0,0]);
            colormap('gray')
            freezeColors()

%             S2=surface(X2,Y2,Z2,'FaceAlpha',0.2);
%             colormap('pink')
%             freezeColors()
% 
%             
%             S3=surface(X3,Y3,Z3,'FaceAlpha',0.2);
%             colormap('gray')
%             freezeColors()
            
            S4{i}=surface(X4,Y4,Z4,'FaceColor',colors(i,:),'FaceAlpha',FaceAlpha(i),...
                'EdgeColor',[0 0 0],'EdgeAlpha',0.25);
%             colormap('cool')
            freezeColors()
            
            i=i+1;

%             grid on
        end

        legend([S4{:}], {'Fy=0.05','Fy=0.06','Fy=0.07', 'Fy=0.08', 'Fy=0.12'},...
            'location', 'best')

%         xlabel(xData)
%         ylabel(yData)
%         zlabel(zData)
%         zlim([0 max(max(z1))*1.1])
            
            xlabel('h [-]')
            ylabel('T [s]')
            zlabel('EART [worker-days]')
            set(gca, 'FontSize', 12) 

        hold on
        
        for i=1:1:(numel(unique(self.(ctData)))-5)

            ind = self.(ctData)==uniqueCt(i);
            xtarget=self.(xData)(ind).*(self.isCANDIDATE(ind,1));
            xtarget(xtarget==0)=NaN;
            ytarget=self.(yData)(ind).*(self.isCANDIDATE(ind,1));
            ytarget(ytarget==0)=NaN;
            ztarget=self.(zData)(ind,1).*(self.isCANDIDATE(ind,1));
            ztarget(ztarget==0)=NaN;
    
            DC=scatter3(xtarget,ytarget,ztarget,...
                    130, 'MarkerFaceColor',[0.8 0 0], ...
                    'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha',0.5, ...
                    'DisplayName', 'Design candidates'); 

            [0.466 0.674 0.188];
        end

        Starg=surface(X1,Y1,self.EALtarget.*ones(length(X1),length(Y1)),...
            'FaceAlpha',0.8,'EdgeAlpha',0,'EdgeColor',[0 0 0],'FaceColor',[0.5 0.5 0.5],...
            'DisplayName', 'Loss Target');

        end
end

function clearVars()

    clearvars -except design fullFit fullFitAcc

end

end  

end

