classdef coloniesDiffusionAndGrowth
    properties %(SetAccess = protected)
        withAutoInduction = 1;
        withAdaptive = 0;
        
        sideLength = 0;
        midPt = 0;
        quarterPt = 0;
        
        DR = .030; % Diffusivity mm^2/min
        DL = .030;
        DC = .0556; %bionumbers
        DN = .045; %http://www.sciencedirect.com/science/article/pii/S0022024800005285
        DI = .045; %assumed to be close to cdag.DN
        
        YC = .65;
        YN=4.1627;
        YI=3895;
        YNi = 2.0011;
        YIi = 4.9892;
        
        kprodR = .001; % Production rate molecules/second
        kprodL = .001; % Production rate molecules/second
        kdecayR = 0;%.0001; %Decay rate molecules/second
        kdecayL = 0;%.0001; %Decay rate molecules/second
        kdecayX = 0;%.0062/60; %hourly rate scaled to minutes
        kautoR = 0.0001;
        kautoL = 0.0001;
        kcrossL = 0.0005;
        
        dx = .1; %units of mm
        dy = .1; %units of mm
        dt = 0;
        tend = 0; % End time (s)
        % cdag.timesToPrint = [1 25 50 75 100 200 300 500 1000];
        timesToPrint = [1 250 500 750 1000 1250 1500 1750 2000];
        %cdag.timesToPrint = cdag.dt*[1 2 3 4 5 6 7 8 9 10];
        cmap = [];
        numTimeSteps = 0;
        frameTimeSteps = [];
        
        X = [];
        XThresh = 10000*1/(.15*1.1*10^9);%.0001
        LIKO = [];
        RIKO = [];
        LRKO = [];
        RRKO = [];
        colonies = [];
        colonyCenterCoords = [];
        R = [];
        L = [];
        initialR = 0;
        initialL = 0;
        Rstore = [];
        Lstore = [];
        GFP = [];
        C = [];
        N = [];
        I = [];
        Ni = [];
        Ii = [];
        mu = [];
        muMax = 0.3341/60; % hourly rate scaled to minutes
        muMaxPrime = .0616/60;
        
        sideOffset=0;
        
        RAll = [];
        LAll = [];
        XAll = [];
        GFPAll = [];
        GFPAll2 = [];
        XAll2 = [];
        CAll2 = [];
        RIKOAll2 = [];
        RAll2 = [];
        
        onPlate = [];
        plateBorderHorz = [];
        plateBorderVert = [];
        pointsBorderToInsideHorz = [];
        pointsBorderToInsideVert = [];
        
        diffKernelHorz = [0 0 0; 1 -2 1; 0 0 0];
        diffKernelVert = [0 1 0; 0 -2 0; 0 1 0];
        
        extendedL1=[];extendedL2=[];
        extendedR1=[];extendedR2=[];
        extendedC1=[];extendedC2=[];
        extendedN1=[];extendedN2=[];
        extendedI1=[];extendedI2=[];
        
        skipsteps=0;
        timeStep=1;
        
        xCum = [];
        cCum = [];
        nCum = [];
        iCum = [];
        lCum = [];
        rCum = [];
    end
    
    methods
        
        function cdag = coloniesDiffusionAndGrowth(aSideLength,aTEnd)
            
            %set time steps
            cdag.tend = aTEnd;
            cdag.dt = cdag.dx^2/max([cdag.DR,cdag.DL,cdag.DC,cdag.DN,cdag.DI])*.2495; %CFL constant
            cdag.numTimeSteps = length(0:cdag.dt:cdag.tend);            
            cdag.frameTimeSteps = 1:ceil(cdag.numTimeSteps/100):cdag.numTimeSteps;
            
            cdag.cmap = jet(length(cdag.timesToPrint));
            
            %set spatial dimensions
            cdag.sideLength = aSideLength;
            cdag.midPt = ceil(cdag.sideLength/2.0);
            cdag.quarterPt = floor(cdag.midPt/2);
            
            %set biomass, genotype, and colonies matrices
            cdag.X = zeros(cdag.sideLength, cdag.sideLength);
            cdag.LIKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.RIKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.LRKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.RRKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.colonies = zeros(cdag.sideLength, cdag.sideLength);
            
            %set producer colony coordinates
            cdag.X(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350), ...
                ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350)) = ...
                100/(.15*1.1*10^11);
            cdag.colonies(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350), ...
                ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350)) = ...
                1;
            cdag.colonyCenterCoords = [ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350) ...
                ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350)];
            
            %set RIKO colony coordinates
            RIKOcoords=[160 25;
                210 25;
                275 65;
                110 75;
                80 100;
                310 100;
                40 120;
                125 120;
                175 120;
                35 135;
                280 135;
                145 150;
                310 150;
                20 160;
                260 170;
                90 190;
                220 200;
                325 210;
                55 220;
                105 225;
                125 250;
                155 275;
                120 290;
                220 295;
                160 300];
            for i=1:length(RIKOcoords)
                cdag.X(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350), ...
                    ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350)) = ...
                    100/(.15*1.1*10^11);
                cdag.RIKO(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350), ...
                    ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350)) = ...
                    1;
                cdag.colonies(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350), ...
                    ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350)) = ...
                    i+1;
                cdag.colonyCenterCoords(end+1,:)=[ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350) ...
                    ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350)];
            end
            
            %set autoinducer, nutrient, and internal nutrient matrices
            cdag.R = zeros(cdag.sideLength, cdag.sideLength);
            cdag.L = zeros(cdag.sideLength, cdag.sideLength);
            cdag.R(cdag.midPt,cdag.midPt) = cdag.initialR;
            cdag.L(cdag.midPt,cdag.midPt) = cdag.initialL;
            cdag.Rstore = zeros(1,cdag.numTimeSteps);
            cdag.Lstore = zeros(1,cdag.numTimeSteps);
            %does this makes sense given nonzero R?
            cdag.GFP = cdag.R.*cdag.X.*cdag.RIKO;
            cdag.C = .5*.0012*1/(cdag.sideLength*cdag.sideLength)*ones(cdag.sideLength, cdag.sideLength);
            cdag.N = .0625*.0012*1/(cdag.sideLength*cdag.sideLength)*ones(cdag.sideLength, cdag.sideLength);
            cdag.I = .00006981*.0012*1/(cdag.sideLength*cdag.sideLength)*ones(cdag.sideLength, cdag.sideLength);
            cdag.Ni = ones(cdag.sideLength, cdag.sideLength).*(cdag.X~=0);
            cdag.Ii = ones(cdag.sideLength, cdag.sideLength).*(cdag.X~=0);
            
            cdag.RAll = zeros(cdag.sideLength, cdag.sideLength, length(cdag.timesToPrint));
            cdag.LAll = zeros(cdag.sideLength, cdag.sideLength, length(cdag.timesToPrint));
            cdag.XAll = zeros(cdag.sideLength, cdag.sideLength, length(cdag.timesToPrint));
            cdag.GFPAll = zeros(cdag.sideLength, cdag.sideLength, length(cdag.timesToPrint));
            cdag.GFPAll2 = zeros(cdag.sideLength, cdag.sideLength, length(cdag.frameTimeSteps));
            cdag.XAll2 = zeros(cdag.sideLength, cdag.sideLength, length(cdag.frameTimeSteps));
            cdag.CAll2 = zeros(cdag.sideLength, cdag.sideLength, length(cdag.frameTimeSteps));
            cdag.RIKOAll2 = zeros(cdag.sideLength, cdag.sideLength, length(cdag.frameTimeSteps));
            cdag.RAll2 = zeros(cdag.sideLength, cdag.sideLength, length(cdag.frameTimeSteps));
            
            %calculate boolean matrix onPlate representing a circular plate
            [coords1 coords2] = meshgrid(1:cdag.sideLength,1:cdag.sideLength);
            cdag.onPlate = (sqrt((coords1-cdag.midPt).^2 + (coords2-cdag.midPt).^2) <= (cdag.midPt-1));
            %calculate plateBorderHorz and Vert as boolean matrices, which
            %represent points not on the plate, but which horizontally or
            %vertically border the plate
            %for actual calculation, take the intersection of all points not
            %on the plate, with horizontally or vertically shifted by one
            %versions of onPlate
            cdag.plateBorderHorz = ~cdag.onPlate & ...
                ([zeros(cdag.sideLength,1) cdag.onPlate(:,1:end-1)] | ...
                [cdag.onPlate(:,2:end) zeros(cdag.sideLength,1)]);
            cdag.plateBorderVert = ~cdag.onPlate & ...
                ([zeros(1,cdag.sideLength); cdag.onPlate(1:end-1,:)] | ...
                [cdag.onPlate(2:end,:); zeros(1,cdag.sideLength)]);
            %We must have a special case whenever onPlate consists of an
            %isolated point in any row or column. For example, consider the
            %following figure, where 1 represents regular onPlate points,
            %and 2 isolated points, 
            % 0 0 2 0 0
            % 0 1 1 1 0
            % 2 1 1 1 2
            % 0 1 1 1 0
            % 0 0 2 0 0
            %In makeExtendedMatrix, we must map the values of isolated
            %points to their horizontal and vertical neighbors. The
            %following code detects isolated points since the corresponding
            %plateBorderHorz or Vert points will be exactly 2 apart, and makes
            %maps pointsBorderToInsideHorz and Vert that take the coordinates of the border
            %points to those of the isolated point
            %It also zeros out the isolated point in the plateBorder
            %matrices, so they will only represent nonisolated points.
            cdag.pointsBorderToInsideHorz = containers.Map('KeyType','double','ValueType','double');
            cdag.pointsBorderToInsideVert = containers.Map('KeyType','double','ValueType','double');
            for i=1:cdag.sideLength
                if(abs(diff(find(cdag.plateBorderHorz(i,:))))==2)
                    jCoords=find(cdag.plateBorderHorz(i,:));
                    cdag.pointsBorderToInsideHorz((jCoords(1)-1)*cdag.sideLength+i)=jCoords(1)*cdag.sideLength+i;
                    cdag.pointsBorderToInsideHorz((jCoords(2)-1)*cdag.sideLength+i)=jCoords(1)*cdag.sideLength+i;
                    cdag.plateBorderHorz(i,:) = zeros(1,cdag.sideLength);
                end
            end
            for i=1:cdag.sideLength
                if(abs(diff(find(cdag.plateBorderVert(:,i))))==2)
                    jCoords=find(cdag.plateBorderVert(:,i));
                    cdag.pointsBorderToInsideVert((i-1)*cdag.sideLength+jCoords(1))=(i-1)*cdag.sideLength+jCoords(1)+1;
                    cdag.pointsBorderToInsideVert((i-1)*cdag.sideLength+jCoords(2))=(i-1)*cdag.sideLength+jCoords(1)+1;
                    cdag.plateBorderVert(:,i) = zeros(cdag.sideLength,1);
                end
            end
            
            %zero all concentration, biomass and genotype matrices not on
            %the plate
            cdag.X(~cdag.onPlate) = 0;
            cdag.R(~cdag.onPlate) = 0;
            cdag.L(~cdag.onPlate) = 0;
            cdag.LIKO(~cdag.onPlate) = 0;
            cdag.RIKO(~cdag.onPlate) = 0;
            cdag.LRKO(~cdag.onPlate) = 0;
            cdag.RRKO(~cdag.onPlate) = 0;
            cdag.C(~cdag.onPlate) = 0;
            cdag.N(~cdag.onPlate) = 0;
            cdag.I(~cdag.onPlate) = 0;
            cdag.Ni(~cdag.onPlate) = 0;
            cdag.Ii(~cdag.onPlate) = 0;
            
            %calculate extended matrices, where the values of the
            %plateBorder points match either their horizontal or vertical
            %onPlate neighbors
            [cdag.extendedL1 cdag.extendedL2] = cdag.makeExtendedMatrix(cdag.L);
            [cdag.extendedR1 cdag.extendedR2] = cdag.makeExtendedMatrix(cdag.R);
            [cdag.extendedC1, cdag.extendedC2] = cdag.makeExtendedMatrix(cdag.C);
            [cdag.extendedN1 cdag.extendedN2] = cdag.makeExtendedMatrix(cdag.N);
            [cdag.extendedI1 cdag.extendedI2] = cdag.makeExtendedMatrix(cdag.I);
        end
        
        %
        %zero values outside the plate to prevent diffusion there
        function diffusionMatrix = makeDiffusionMatrix(cdag,Dmatrix,extendedMatrixHorz,extendedMatrixVert)
            diffusionMatrix = Dmatrix*filter2(cdag.diffKernelHorz,extendedMatrixHorz,'valid')/(cdag.dx)^2 + ...
                Dmatrix*filter2(cdag.diffKernelVert,extendedMatrixVert,'valid')/(cdag.dx)^2;
            diffusionMatrix(~cdag.onPlate) = 0;
        end
        
        function [extendedMatrixHorz extendedMatrixVert] = makeExtendedMatrix(cdag,matrix)
            %Calculate the horizontal and vertical neighbors of the
            %plateBorder points, by taking the intersection of onPlate and
            %horizontally and vertically shifted versions of the plateBorders 
            matrixHorz=matrix;
            matrixHorz(cdag.plateBorderHorz) = matrix(cdag.onPlate & ...
                ([zeros(cdag.sideLength,1) cdag.plateBorderHorz(:,1:end-1)] | ...
                [cdag.plateBorderHorz(:,2:end) zeros(cdag.sideLength,1)]));
            matrixVert=matrix;
            matrixVert(cdag.plateBorderVert) = matrix(cdag.onPlate & ...
                ([zeros(1,cdag.sideLength); cdag.plateBorderVert(1:end-1,:)] | ...
                [cdag.plateBorderVert(2:end,:); zeros(1,cdag.sideLength)]));
            
            %map the values of the isolated points to their horizontal and
            %vertical neighbors
            pointsBorderHorz = keys(cdag.pointsBorderToInsideHorz);
            for i=1:length(pointsBorderHorz)
                matrixHorz(pointsBorderHorz{i}) = matrixHorz(cdag.pointsBorderToInsideHorz(pointsBorderHorz{i}));
            end
            pointsBorderVert = keys(cdag.pointsBorderToInsideVert);
            for i=1:length(pointsBorderVert)
                matrixVert(pointsBorderVert{i}) = matrixVert(cdag.pointsBorderToInsideVert(pointsBorderVert{i}));
            end
            
            %extend the horizontal and vertical borders by one at each end,  
            extendedMatrixHorz = [matrixHorz(1,1) matrixHorz(1,:) matrixHorz(1,end); ...
                matrixHorz(:,1) matrixHorz matrixHorz(:,end); ...
                matrixHorz(end,1) matrixHorz(end,:) matrixHorz(end,end)];
            extendedMatrixVert = [matrixVert(1,1) matrixVert(1,:) matrixVert(1,end); ...
                matrixVert(:,1) matrixVert matrixVert(:,end); ...
                matrixVert(end,1) matrixVert(end,:) matrixVert(end,end)];
        end
        
        function [cdag newOverflowCoordVert newOverflowCoordHorz newOverflowXVal newOverflowLIKOVal newOverflowRIKOVal ...
                newOverflowLRKOVal newOverflowRRKOVal newOverflowIiVal newOverflowNiVal newOverflowColoniesVal] ...
                = overflowCell(varargin)
            
            cdag = varargin{1};

            overflowCoordVert = varargin{2};
            overflowCoordHorz = varargin{3};
            colony = varargin{4};
            
            %check if this is a newly dividing cell, or we have overflow
            %biomass coming in from an adjacent cell
            first=1;
            if(length(varargin)>4)
                overflowXVal = varargin{5};
                overflowLIKOVal = varargin{6};
                overflowRIKOVal = varargin{7};
                overflowLRKOVal = varargin{8};
                overflowRRKOVal = varargin{9};
                overflowIiVal = varargin{10};
                overflowNiVal = varargin{11};
                overflowColoniesVal = varargin{12};
                first=0;
            end
            if(colony~=0)
                centerCoord1 = cdag.colonyCenterCoords(colony,1);
                centerCoord2 = cdag.colonyCenterCoords(colony,2);
            end
            
            newOverflowLIKOVal = cdag.LIKO(overflowCoordVert,overflowCoordHorz);
            newOverflowRIKOVal = cdag.RIKO(overflowCoordVert,overflowCoordHorz);
            newOverflowLRKOVal = cdag.LRKO(overflowCoordVert,overflowCoordHorz);
            newOverflowRRKOVal = cdag.RRKO(overflowCoordVert,overflowCoordHorz);
            newOverflowColoniesVal = cdag.colonies(overflowCoordVert,overflowCoordHorz);
            if(first)
                newOverflowIiVal = cdag.Ii(overflowCoordVert,overflowCoordHorz);
                newOverflowNiVal = cdag.Ni(overflowCoordVert,overflowCoordHorz);
                newOverflowXVal = cdag.X(overflowCoordVert,overflowCoordHorz)/2;
            else
                newOverflowIiVal = cdag.Ii(overflowCoordVert,overflowCoordHorz);
                newOverflowNiVal = cdag.Ni(overflowCoordVert,overflowCoordHorz);
                newOverflowXVal = cdag.X(overflowCoordVert,overflowCoordHorz);
            end
            
            if(1)
                
                %calculate nearest border point, convert to horizontal and
                %vertical coordinates
                [XDist XNearestBorder] = bwdist(cdag.X==0 & cdag.onPlate,'chessboard');
                nearestBorderCoordVert = ...
                    mod(XNearestBorder(overflowCoordVert,overflowCoordHorz)-1,cdag.sideLength)+1;
                nearestBorderCoordHorz = ...
                    floor((double(XNearestBorder(overflowCoordVert,overflowCoordHorz)) -1 )/ (cdag.sideLength*1.0))+1;
                randCoordChanges=[];
                
                %add to randCoordChanges according to whether we have
                %border points in any direction
                if(nearestBorderCoordVert < overflowCoordVert)
                    randCoordChanges(end+1,:) = [-1 0];
                end
                if(nearestBorderCoordVert > overflowCoordVert)
                    randCoordChanges(end+1,:) = [1 0];
                end
                if(nearestBorderCoordHorz < overflowCoordHorz)
                    randCoordChanges(end+1,:) = [0 -1];
                end
                if(nearestBorderCoordHorz > overflowCoordHorz)
                    randCoordChanges(end+1,:) = [0 1];
                end
                if(nearestBorderCoordVert < overflowCoordVert && nearestBorderCoordHorz < overflowCoordHorz)
                    randCoordChanges(end+1,:) = [-1 -1];
                end
                if(nearestBorderCoordVert < overflowCoordVert && nearestBorderCoordHorz > overflowCoordHorz)
                    randCoordChanges(end+1,:) = [-1 1];
                end
                if(nearestBorderCoordVert > overflowCoordVert && nearestBorderCoordHorz < overflowCoordHorz)
                    randCoordChanges(end+1,:) = [1 -1];
                end
                if(nearestBorderCoordVert > overflowCoordVert && nearestBorderCoordHorz > overflowCoordHorz)
                    randCoordChanges(end+1,:) = [1 1];
                end
                
                %if randCoordChanges is not empty, select one of the randCoordChanges 
                %entries, calculate newOverflowCoords with cyclic boundary conditions
                if(length(randCoordChanges)==0)
                    newOverflowCoordVert=overflowCoordVert;
                    newOverflowCoordHorz=overflowCoordHorz;
                else
                    numRand = randi(size(randCoordChanges,1),1);
                    newOverflowCoordVert=overflowCoordVert+randCoordChanges(numRand,1);
                    if(newOverflowCoordVert>cdag.sideLength)
                        newOverflowCoordVert=1;
                    end
                    if(newOverflowCoordVert<1)
                        newOverflowCoordVert=cdag.sideLength;
                    end
                    newOverflowCoordHorz=overflowCoordHorz+randCoordChanges(numRand,2);
                    if(newOverflowCoordHorz>cdag.sideLength)
                        newOverflowCoordHorz=1;
                    end
                    if(newOverflowCoordHorz<1)
                        newOverflowCoordHorz=cdag.sideLength;
                    end
                end
            end
            
            if(first)
                cdag.Ii(overflowCoordVert,overflowCoordHorz) = cdag.Ii(overflowCoordVert,overflowCoordHorz);
                cdag.Ni(overflowCoordVert,overflowCoordHorz) = cdag.Ni(overflowCoordVert,overflowCoordHorz);
                cdag.X(overflowCoordVert,overflowCoordHorz) = cdag.X(overflowCoordVert,overflowCoordHorz)/2;
            else
                cdag.LIKO(overflowCoordVert,overflowCoordHorz) = overflowLIKOVal;
                cdag.RIKO(overflowCoordVert,overflowCoordHorz) = overflowRIKOVal;
                cdag.LRKO(overflowCoordVert,overflowCoordHorz) = overflowLRKOVal;
                cdag.RRKO(overflowCoordVert,overflowCoordHorz) = overflowRRKOVal;
                cdag.colonies(overflowCoordVert,overflowCoordHorz) = overflowColoniesVal;
                cdag.Ii(overflowCoordVert,overflowCoordHorz) = overflowIiVal;
                cdag.Ni(overflowCoordVert,overflowCoordHorz) = overflowNiVal;
                cdag.X(overflowCoordVert,overflowCoordHorz) = overflowXVal;
            end
        end
        
        function main(cdag)
            tic
            while(cdag.timeStep <= cdag.numTimeSteps)
                cdag = cdag.runOneTimeStep();
            end
            toc
            cdag = cdag.makeGraphics();
        end
        
        function cdag = runOneTimeStep(cdag)
            disp(num2str(cdag.timeStep*cdag.dt))
            
            cdag.Lstore(cdag.timeStep) = cdag.L(cdag.midPt,cdag.midPt);
            cdag.Rstore(cdag.timeStep) = cdag.R(cdag.midPt,cdag.midPt);

            %L dynamics
            diffusionTermL = cdag.makeDiffusionMatrix(cdag.DL,cdag.extendedL1,cdag.extendedL2);
            cdag.L = cdag.L + cdag.dt*diffusionTermL;
            productionTermL = cdag.kprodL*(cdag.X.*(~cdag.LIKO));
            if(cdag.withAutoInduction)
                productionTermL = productionTermL + ...
                    cdag.X.*(~cdag.LIKO).*(~cdag.LRKO).* ...
                    cdag.kprodL.*cdag.L./(cdag.kautoL+cdag.L);
            end
            cdag.L = cdag.L + cdag.dt*productionTermL;
            decayTermL = cdag.kdecayL*cdag.L;
            cdag.L = cdag.L + cdag.dt*decayTermL;
            %changeTermL = diffusionTermL + productionTermL - decayTermL;
            
            %R dynamics
            diffusionTermR = cdag.makeDiffusionMatrix(cdag.DR,cdag.extendedR1,cdag.extendedR2);
            cdag.R = cdag.R + cdag.dt*diffusionTermR;
            productionTermR = cdag.kprodR*(cdag.X.*(~cdag.RIKO));
            if(cdag.withAutoInduction)
                productionTermR = productionTermR + ...
                    (cdag.X.*(~cdag.RIKO).*(~cdag.RRKO).*(~cdag.LRKO).* ...
                    cdag.kprodR.*cdag.R./(cdag.kautoR+cdag.R).* (cdag.L./(cdag.kcrossL + cdag.L)));
            end
            cdag.R = cdag.R + cdag.dt*productionTermR;
            decayTermR = cdag.kdecayR*cdag.R;
            cdag.R = cdag.R + cdag.dt*decayTermR;
            %changeTermR = diffusionTermR + productionTermR - decayTermR;
            
            %nutrient diffusion
            diffusionTermC = cdag.makeDiffusionMatrix(cdag.DC,cdag.extendedC1,cdag.extendedC2);
            diffusionTermN = cdag.makeDiffusionMatrix(cdag.DN,cdag.extendedN1,cdag.extendedN2);
            diffusionTermI = cdag.makeDiffusionMatrix(cdag.DI,cdag.extendedI1,cdag.extendedI2);
            changeTermC = diffusionTermC;
            changeTermN = diffusionTermN;
            changeTermI = diffusionTermI;
            cdag.C = max(0,cdag.C + cdag.dt*changeTermC);
            cdag.N = max(0,cdag.N + cdag.dt*changeTermN);
            cdag.I = max(0,cdag.I + cdag.dt*changeTermI);
            
            %growth rate
            cdag.mu = cdag.muMax*(cdag.C>0 & cdag.N>0).*cdag.Ii + cdag.muMaxPrime*(cdag.N<=0 & cdag.C>0).*cdag.Ni;
            growthTermX = cdag.mu.*cdag.X.*(cdag.C>0) + (cdag.mu-cdag.kdecayX).*cdag.X.*(cdag.C==0);
            changeTermX = growthTermX;
            
            %nutrient consumptions
            consumptionTermC = -1/cdag.YC*cdag.mu.*cdag.X.*(cdag.C>0);
            consumptionTermN = -1/cdag.YN*cdag.mu.*cdag.X.*(cdag.N>0);
            consumptionTermI = -1/cdag.YI*cdag.mu.*cdag.X.*(cdag.I>0);
            consumptionTermNi = -(1/cdag.YNi+cdag.Ni).*cdag.mu.*(cdag.N==0 & cdag.Ni>0);
            consumptionTermIi = -(1/cdag.YIi+cdag.Ii).*cdag.mu.*(cdag.I==0 & cdag.Ii>0);
            changeTermC = consumptionTermC;
            changeTermN = consumptionTermN;
            changeTermI = consumptionTermI;
            changeTermNi = consumptionTermNi;
            changeTermIi = consumptionTermIi;
            cdag.C = max(0,cdag.C + cdag.dt*changeTermC);
            cdag.N = max(0,cdag.N + cdag.dt*changeTermN);
            cdag.I = max(0,cdag.I + cdag.dt*changeTermI);
            cdag.Ni = max(0,cdag.Ni + cdag.dt*changeTermNi);
            cdag.Ii = max(0,cdag.Ii + cdag.dt*changeTermIi);
            cdag.X = max(0,cdag.X + cdag.dt*changeTermX);
            
            %zero nutrient concentrations that are below the minimum
            %consumption of one cell, divided by the points on the grid
            if(min(min(consumptionTermC))~=0)
                cdag.C(cdag.C<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermC(consumptionTermC~=0))))=0;
            end
            if(min(min(consumptionTermN))~=0)
                cdag.N(cdag.N<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermN(consumptionTermN~=0))))=0;
            end
            if(min(min(consumptionTermI))~=0)
                cdag.I(cdag.I<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermI(consumptionTermI~=0))))=0;
            end
            
            %find any cell that is above X threshold, divide it in half and
            %move biomass somewhere else. If there is already a cell at
            %that biomass, keep moving, always in the direction of colony
            %border, till you hit an empty grid position
            [overflowCoordVert overflowCoordHorz] = find(cdag.X>cdag.XThresh,1,'first');
            while(~isempty(overflowCoordVert))
                [cdag overflowCoordVert,overflowCoordHorz,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
                    overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal] = ...
                    cdag.overflowCell(overflowCoordVert,overflowCoordHorz, ...
                    cdag.colonies(overflowCoordVert,overflowCoordHorz));
                while( overflowXVal > 0 )
                    [cdag overflowCoordVert,overflowCoordHorz,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
                        overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal] = ...
                        cdag.overflowCell(overflowCoordVert,overflowCoordHorz, ...
                        cdag.colonies(overflowCoordVert,overflowCoordHorz), ...
                        overflowXVal,overflowLIKOVal,overflowRIKOVal,overflowLRKOVal, ...
                        overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal);
                end
                [overflowCoordVert overflowCoordHorz] = find(cdag.X>cdag.XThresh,1,'first');
            end
            
            %calculate cumulative values
            cdag.xCum(cdag.timeStep) = sum(sum(cdag.X));
            cdag.cCum(cdag.timeStep) = sum(sum(cdag.C));
            cdag.nCum(cdag.timeStep) = sum(sum(cdag.N));
            cdag.iCum(cdag.timeStep) = sum(sum(cdag.I));
            cdag.rCum(cdag.timeStep) = sum(sum(cdag.R));
            cdag.lCum(cdag.timeStep) = sum(sum(cdag.L));
            
            %store R, L, X, and GFP for either printing or movies
            matchesTimesToPrint = (uint64(cdag.timesToPrint/cdag.dt) == cdag.timeStep);
            if sum(matchesTimesToPrint)~=0
                cdag.RAll(:,:,matchesTimesToPrint) = cdag.R;
                cdag.LAll(:,:,matchesTimesToPrint) = cdag.L;
                cdag.XAll(:,:,matchesTimesToPrint) = cdag.X;
                cdag.GFPAll(:,:,matchesTimesToPrint) = cdag.GFP;
            end
            
            if sum(cdag.frameTimeSteps == cdag.timeStep)~=0
                cdag.GFPAll2(:,:,cdag.frameTimeSteps == cdag.timeStep) = cdag.GFP;
                cdag.CAll2(:,:,cdag.frameTimeSteps == cdag.timeStep) = cdag.C;
                cdag.XAll2(:,:,cdag.frameTimeSteps == cdag.timeStep) = cdag.X;
                cdag.RIKOAll2(:,:,cdag.frameTimeSteps == cdag.timeStep) = cdag.RIKO;
                cdag.RAll2(:,:,cdag.frameTimeSteps == cdag.timeStep) = cdag.R;
            end
            
            %increment L, R, and GFP
            if(cdag.withAdaptive)
                [Lnext newTimeStepL] = adaptiveTimeStep(cdag.L,cdag.dt,changeTermL,cdag.timeStep);
                [Rnext newTimeStepR] = adaptiveTimeStep(cdag.R,cdag.dt,changeTermR,cdag.timeStep);
                Lnext = cdag.L + min(newTimeStepL, newTimeStepR) * (changeTermL);
                Rnext = cdag.R + min(newTimeStepL, newTimeStepR) * (changeTermR);
            else
                muTemp1 = zeros(cdag.sideLength,cdag.sideLength);
                muTemp2 = zeros(cdag.sideLength,cdag.sideLength);
                muTemp1(cdag.mu~=0) = cdag.muMaxPrime./cdag.mu(cdag.mu~=0)-1;
                muTemp1(muTemp1<0) = 0;
                muTemp2(cdag.mu~=0) = cdag.muMax./cdag.mu(cdag.mu~=0)-1;
                muTemp2(muTemp2<0) = 0;
                changeTermGFP = 18542*cdag.X + ...
                    ( 7454.1*(muTemp1)*(1./(1 + (.1667*muTemp2).^5.8725)) ).*(cdag.N==0).*(cdag.X~=0) + ...
                    ( 8593.3*(muTemp2)*(1./(1 + (.3293*muTemp2).^4.9892)) ).*(cdag.I==0).*(cdag.X~=0);
                cdag.GFP = cdag.GFP + cdag.dt*changeTermGFP.*cdag.RIKO.*cdag.R;
                cdag.timeStep = cdag.timeStep + 1;
            end            
            
            %make new extended matrices
            [cdag.extendedL1 cdag.extendedL2] = cdag.makeExtendedMatrix(cdag.L);
            [cdag.extendedR1 cdag.extendedR2] = cdag.makeExtendedMatrix(cdag.R);
            [cdag.extendedC1 cdag.extendedC2] = cdag.makeExtendedMatrix(cdag.C);
            [cdag.extendedN1 cdag.extendedN2] = cdag.makeExtendedMatrix(cdag.N);
            [cdag.extendedI1 cdag.extendedI2] = cdag.makeExtendedMatrix(cdag.I);
        end
        
        function cdag = makeGraphics(cdag)
            %make three figures for cdag.L, cdag.R, and combined images
            fig1 = figure(1);
            hold on
            fig2 = figure(2);
            hold on
            fig3 = figure(3);
            hold on
            fig4 = figure(4);
            hold on
            fig11 = figure(11);
            hold on
            
            cdag.timeStep = 1;
            
            figArray = [fig1, fig2, fig3, fig4, fig11];
            
            matrixToPlot = zeros(cdag.sideLength,cdag.sideLength,3);
            for i = 1:5
                set(0, 'CurrentFigure', figArray(i))
                
                %find min and max of cdag.L, cdag.R, or both across all time points to plot
                if(i==1)
                    MINVAL=min(min(min(cdag.LAll)));
                    MAXVAL=max(max(max(cdag.LAll)));
                elseif(i==2)
                    MINVAL=min(min(min(cdag.RAll)));
                    MAXVAL=max(max(max(cdag.RAll)));
                elseif(i==3)
                    MINVAL=min(min(min(min(cdag.RAll))),min(min(min(cdag.LAll))));
                    MAXVAL=max(max(max(max(cdag.RAll))),max(max(max(cdag.LAll))));
                elseif(i==4)
                    MINVAL=min(min(min(cdag.XAll)));
                    MAXVAL=max(max(max(cdag.XAll)));
                else
                    MINVAL=min(min(min(cdag.GFPAll)));
                    MAXVAL=max(max(max(cdag.GFPAll)));
                end
                
                for j=1:length(cdag.timesToPrint)
                    cdag.R = cdag.RAll(:,:,j);
                    cdag.L = cdag.LAll(:,:,j);
                    cdag.X = cdag.XAll(:,:,j);
                    cdag.GFP = cdag.GFPAll(:,:,j);
                    
                    %code cdag.L as blue, cdag.R as green
                    if i == 1
                        matrixToPlot(:, :, 1) = uint16(0);
                        matrixToPlot(:, :, 2) = uint16(0);
                        matrixToPlot(:, :, 3) = cdag.L;
                    elseif i==2
                        matrixToPlot(:, :, 1) = uint16(0);
                        matrixToPlot(:, :, 2) = cdag.R;
                        matrixToPlot(:, :, 3) = uint16(0);
                    elseif i==3
                        matrixToPlot(:, :, 1) = uint16(0);
                        matrixToPlot(:, :, 2) = cdag.R;
                        matrixToPlot(:, :, 3) = cdag.L;
                    elseif i==4
                        matrixToPlot(:, :, 1) = cdag.X;
                        matrixToPlot(:, :, 2) = uint16(0);
                        matrixToPlot(:, :, 3) = uint16(0);
                    else
                        matrixToPlot(:, :, 1) = uint16(0);
                        matrixToPlot(:, :, 2) = cdag.GFP;
                        matrixToPlot(:, :, 3) = uint16(0);
                    end
                    
                    %scale images by appropriate maximum across all time points to
                    %plot
                    matrixToPlotScale = uint16( (matrixToPlot - MINVAL) / (MAXVAL - MINVAL) * (2^16-1) );
                    if(i==5)
                        matrixToPlotScale(:,:,1) = cdag.X.*(~cdag.RIKO)*(2^16-1);
                    end
                    
                    %plot horizontal cross-section for cdag.L and cdag.R images only
                    subplot(ceil((length(cdag.timesToPrint)+1)/4),4,1)
                    set(gca,'XLim',[1 cdag.sideLength])
                    if i==1
                        plot(1:cdag.sideLength, matrixToPlot(cdag.midPt,:,3), 'Color', ...
                            cdag.cmap(j,:));
                        xlabel('x')
                        ylabel('cdag.L')
                    elseif i==2
                        plot(1:cdag.sideLength, matrixToPlot(cdag.midPt,:,2), 'Color', ...
                            cdag.cmap(j,:));
                        xlabel('x')
                        ylabel('cdag.R')
                    end
                    hold on
                    
                    %plot 2d images of cdag.L, cdag.R, or both
                    if i==1 || i==2
                        subplot(ceil((length(cdag.timesToPrint)+1)/4), ...
                            4, 1 + j)
                    else
                        subplot(ceil((length(cdag.timesToPrint))/4), ...
                            4, j)
                    end
                    imshow(matrixToPlotScale)
                    title(['t = ' num2str(cdag.timesToPrint(j))])
                    axis([1 cdag.sideLength 1 cdag.sideLength])
                    hold on
                end
            end
            
            figure(10)
            hold on
            subplot(3,4,1)
            plot(1:cdag.numTimeSteps, cdag.xCum, 'k')
            subplot(3,4,2)
            plot(1:cdag.numTimeSteps, cdag.cCum, 'm')
            subplot(3,4,3)
            plot(1:cdag.numTimeSteps, cdag.nCum, 'c')
            subplot(3,4,4)
            plot(1:cdag.numTimeSteps, cdag.iCum, 'y')
            subplot(3,4,5)
            plot(1:cdag.numTimeSteps, cdag.rCum, 'g')
            subplot(3,4,6)
            imshow(cdag.C)
            subplot(3,4,7)
            imshow(cdag.N)
            subplot(3,4,8)
            imshow(cdag.I)
            subplot(3,4,9)
            plot(1:cdag.numTimeSteps, cdag.lCum, 'b')
            
            figure(5)
            surf(cdag.C)
            figure(6)
            surf(cdag.N)
            figure(7)
            surf(cdag.I)
            figure(8)
            surf(cdag.Ni)
            figure(9)
            surf(cdag.Ii)
            
            if(~exist('frames','dir'))
                eval('mkdir frames');
            end
            
            MINVALGFP = min(min(min(cdag.GFPAll2)));
            MAXVALGFP = max(max(max(cdag.GFPAll2)));
            MINVALX = min(min(min(cdag.XAll2)));
            MAXVALX = max(max(max(cdag.XAll2)));
            MINVALR = min(min(min(cdag.RAll2)));
            MAXVALR = max(max(max(cdag.RAll2)));
            for i=1:length(cdag.frameTimeSteps)
                cdag.GFP = cdag.GFPAll2(:,:,i);
                %disp(max(max(cdag.GFPAll2(:,:,i))))
                cdag.X = cdag.XAll2(:,:,i);
                cdag.RIKO = cdag.RIKOAll2(:,:,i);
                cdag.R = cdag.RAll2(:,:,i);
                %cdag.GFP = cdag.GFP + cdag.RIKO*.1*MAXVAL1;
                matrixToPlot(:, :, 1) = cdag.X.*(~cdag.RIKO);
                matrixToPlot(:, :, 2) = cdag.GFP;
                matrixToPlot(:, :, 3) = cdag.R.*(cdag.X==0);
                matrixToPlotScale = uint16(zeros(size(matrixToPlot,1),size(matrixToPlot,2),size(matrixToPlot,3)));
                matrixToPlotScale(:, :, 3) = uint16( (matrixToPlot(:, :, 3) - MINVALR) / (MAXVALR - MINVALR) * (2^16-1) );
                matrixToPlotScale(:, :, 2) = uint16( (matrixToPlot(:, :, 2) - MINVALGFP) / (MAXVALGFP - MINVALGFP) * (2^16-1) );
                %add baseline GFP value corresponding to 2000 on the RGB
                %scale
                matrixToPlotScale(:, :, 2) = min( uint16(matrixToPlotScale(:, :, 2) + uint16(cdag.RIKO*2000)), uint16(2^16-1));
                matrixToPlotScale(:, :, 1) = uint16( (matrixToPlot(:, :, 1) - MINVALX) / (MAXVALX - MINVALX) * (2^16-1) );
                figure(15)
                imshow(matrixToPlotScale ...
                    (cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,:))
                drawnow
                print( 15, '-djpeg', sprintf('frames/frmGFP%d.jpg', i))
            end
            
            disp('making ffmpeg movie')
            system('"C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmGFP%d.jpg frames/GFP.mov')
            disp('done with ffmpeg movie')
            
            MINVALC=min(min(min(cdag.CAll2)));
            MAXVALC=max(max(max(cdag.CAll2)));
            MINVALX=min(min(min(cdag.XAll2)));
            MAXVALX=max(max(max(cdag.XAll2)));
            for i=1:length(cdag.frameTimeSteps)
                cdag.C = cdag.CAll2(:,:,i);
                cdag.X = cdag.XAll2(:,:,i);
                matrixToPlot(:, :, 1) = cdag.X;
                matrixToPlot(:, :, 2) = uint16(0);
                matrixToPlot(:, :, 3) = cdag.C;
                matrixToPlotScale = uint16(zeros(size(matrixToPlot,1),size(matrixToPlot,2),size(matrixToPlot,3)));
                matrixToPlotScale(:,:,3) = uint16( (matrixToPlot(:,:,3) - MINVALC) / (MAXVALC - MINVALC) * (2^16-1) );
                matrixToPlotScale(:,:,1) = uint16( (matrixToPlot(:,:,1) - MINVALX) / (MAXVALX - MINVALX) * (2^16-1) );
                figure(15)
                imshow(matrixToPlotScale ...
                    (cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,:))
                drawnow
                print( 15, '-djpeg', sprintf('frames/frmX%d.jpg', i))
            end
            
            disp('making ffmpeg movie')
            system('"C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmX%d.jpg frames/X.mov')
            disp('done with ffmpeg movie')
            %end
        end
    end
end