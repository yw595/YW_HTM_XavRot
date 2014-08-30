classdef coloniesDiffusionAndGrowth
    properties %(SetAccess = protected)
        withAutoInduction = 1;
        withAdaptive = 0;
        
        sideLength = 200;
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
        tend = 2400; % End time (s)
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
        horzPlateEdge = [];
        vertPlateEdge = [];
        pointsBorderToInside1 = [];
        pointsBorderToInside2 = [];
        
        diffKernel1 = [0 0 0; 1 -2 1; 0 0 0];
        diffKernel2 = [0 1 0; 0 -2 0; 0 1 0];
        
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
        
        function cdag = coloniesDiffusionAndGrowth()
            
            cdag.dt = cdag.dx^2/max([cdag.DR,cdag.DL,cdag.DC,cdag.DN,cdag.DI])*.2495; %CFL constant
            
            cdag.cmap = jet(length(cdag.timesToPrint));
            cdag.numTimeSteps = length(0:cdag.dt:cdag.tend);
            cdag.frameTimeSteps = 1:ceil(cdag.numTimeSteps/100):cdag.numTimeSteps;
            
            cdag.midPt = ceil(cdag.sideLength/2.0);
            cdag.quarterPt = floor(cdag.midPt/2);
            
            cdag.X = zeros(cdag.sideLength, cdag.sideLength);
            cdag.LIKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.RIKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.LRKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.RRKO = zeros(cdag.sideLength, cdag.sideLength);
            cdag.colonies = zeros(cdag.sideLength, cdag.sideLength);
            
            cdag.X(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350),ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350))=100/(.15*1.1*10^11);
            cdag.colonies(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350),ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350))=1;
            cdag.colonyCenterCoords = [ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*100/350) ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*50/350)];
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
                cdag.X(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350),ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350))=100/(.15*1.1*10^11);
                cdag.RIKO(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350),ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350))=1;
                cdag.colonies(ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350),ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350))=i+1;
                cdag.colonyCenterCoords(end+1,:)=[ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,1)/350) ceil(cdag.sideOffset+(cdag.sideLength-2*cdag.sideOffset)*RIKOcoords(i,2)/350)];
            end
            
            cdag.R = zeros(cdag.sideLength, cdag.sideLength);
            cdag.L = zeros(cdag.sideLength, cdag.sideLength);
            
            cdag.Rstore = zeros(1,cdag.numTimeSteps);
            cdag.Lstore = zeros(1,cdag.numTimeSteps);
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
            
            cdag.R(cdag.midPt,cdag.midPt) = cdag.initialR;
            cdag.L(cdag.midPt,cdag.midPt) = cdag.initialL;
            
            [coords1 coords2] = meshgrid(1:cdag.sideLength,1:cdag.sideLength);
            cdag.onPlate = (sqrt((coords1-cdag.midPt).^2 + (coords2-cdag.midPt).^2) <= (cdag.midPt-1));
            cdag.horzPlateEdge = ~cdag.onPlate & ([zeros(cdag.sideLength,1) cdag.onPlate(:,1:end-1)] | [cdag.onPlate(:,2:end) zeros(cdag.sideLength,1)]);
            cdag.vertPlateEdge = ~cdag.onPlate & ([zeros(1,cdag.sideLength); cdag.onPlate(1:end-1,:)] | [cdag.onPlate(2:end,:); zeros(1,cdag.sideLength)]);
            cdag.pointsBorderToInside1 = containers.Map('KeyType','double','ValueType','double');
            cdag.pointsBorderToInside2 = containers.Map('KeyType','double','ValueType','double');
            for i=1:cdag.sideLength
                if(abs(diff(find(cdag.horzPlateEdge(i,:))))==2)
                    jCoords=find(cdag.horzPlateEdge(i,:));
                    cdag.pointsBorderToInside1((jCoords(1)-1)*cdag.sideLength+i)=jCoords(1)*cdag.sideLength+i;
                    cdag.pointsBorderToInside1((jCoords(2)-1)*cdag.sideLength+i)=jCoords(1)*cdag.sideLength+i;
                    %else
                    cdag.horzPlateEdge(i,:) = zeros(1,cdag.sideLength);
                end
            end
            for i=1:cdag.sideLength
                if(abs(diff(find(cdag.vertPlateEdge(:,i))))==2)
                    jCoords=find(cdag.vertPlateEdge(:,i));
                    cdag.pointsBorderToInside2((i-1)*cdag.sideLength+jCoords(1))=(i-1)*cdag.sideLength+jCoords(1)+1;
                    cdag.pointsBorderToInside2((i-1)*cdag.sideLength+jCoords(2))=(i-1)*cdag.sideLength+jCoords(1)+1;
                    %else
                    cdag.vertPlateEdge(:,i) = zeros(cdag.sideLength,1);
                end
            end
            
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
            
            [cdag.extendedL1 cdag.extendedL2] = cdag.makeExtendedMatrix(cdag.L);
            [cdag.extendedR1 cdag.extendedR2] = cdag.makeExtendedMatrix(cdag.R);
            [cdag.extendedC1 cdag.extendedC2] = cdag.makeExtendedMatrix(cdag.C);
            [cdag.extendedN1 cdag.extendedN2] = cdag.makeExtendedMatrix(cdag.N);
            [cdag.extendedI1 cdag.extendedI2] = cdag.makeExtendedMatrix(cdag.I);
        end
        
        function diffusionMatrix = makeDiffusionMatrix(cdag,matrix,Dmatrix,extendedMatrix1,extendedMatrix2)
            diffusionMatrix = Dmatrix*filter2(cdag.diffKernel1,extendedMatrix1,'valid')/(cdag.dx)^2 + ...
                Dmatrix*filter2(cdag.diffKernel2,extendedMatrix2,'valid')/(cdag.dx)^2;
            diffusionMatrix(~cdag.onPlate) = 0;
            if(0)
                cdag.midPt=size(matrix,1)/2;
                diffusionMatrix = zeros(size(matrix,1),size(matrix,2));
                diffusionMatrix(2:end-1,2:end-1) = Dmatrix*(matrix(2:end-1,1:end-2) + matrix(2:end-1,3:end) - ...
                    4*matrix(2:end-1,2:end-1) + matrix(1:end-2,2:end-1) + matrix(3:end,2:end-1))/(cdag.dx)^2;
                diffusionMatrix(1,2:end-1) = Dmatrix*(matrix(1,1:end-2) + matrix(1,3:end) + matrix(2,2:end-1) - 3*matrix(1,2:end-1))/(cdag.dx)^2;
                diffusionMatrix(end,2:end-1) = Dmatrix*(matrix(end,1:end-2) + matrix(end,3:end) + matrix(end-1,2:end-1) - 3*matrix(end,2:end-1))/(cdag.dx)^2;
                diffusionMatrix(2:end-1,1) = Dmatrix*(matrix(1:end-2,1) + matrix(3:end,1) + matrix(2:end-1,2) - 3*matrix(2:end-1,1))/(cdag.dx)^2;
                diffusionMatrix(2:end-1,end) = Dmatrix*(matrix(1:end-2,end) + matrix(3:end,end) + matrix(2:end-1,end-1) - 3*matrix(2:end-1,end))/(cdag.dx)^2;
                diffusionMatrix(1,1) = Dmatrix*(matrix(2,1) + matrix(1,2) - 2*matrix(1,1))/(cdag.dx)^2;
                diffusionMatrix(1,end) = Dmatrix*(matrix(2,end) + matrix(1,end-1) - 2*matrix(1,end))/(cdag.dx)^2;
                diffusionMatrix(end,1) = Dmatrix*(matrix(end-1,1) + matrix(end,2) - 2*matrix(end,1))/(cdag.dx)^2;
                diffusionMatrix(end,end) = Dmatrix*(matrix(end-1,end) + matrix(end,end-1) - 2*matrix(end,end))/(cdag.dx)^2;
            end
        end
        
        function [extendedMatrix1 extendedMatrix2] = makeExtendedMatrix(cdag,matrix)
            matrix1=matrix;
            matrix1(cdag.horzPlateEdge) = matrix(cdag.onPlate & ([zeros(cdag.sideLength,1) cdag.horzPlateEdge(:,1:end-1)] | ...
                [cdag.horzPlateEdge(:,2:end) zeros(cdag.sideLength,1)]));
            matrix2=matrix;
            matrix2(cdag.vertPlateEdge) = matrix(cdag.onPlate & ([zeros(1,cdag.sideLength); cdag.vertPlateEdge(1:end-1,:)] | ...
                [cdag.vertPlateEdge(2:end,:); zeros(1,cdag.sideLength)]));
            
            pointsBorder1 = keys(cdag.pointsBorderToInside1);
            for i=1:length(pointsBorder1)
                matrix1(pointsBorder1{i}) = matrix1(cdag.pointsBorderToInside1(pointsBorder1{i}));
            end
            pointsBorder2 = keys(cdag.pointsBorderToInside2);
            for i=1:length(pointsBorder2)
                matrix2(pointsBorder2{i}) = matrix2(cdag.pointsBorderToInside2(pointsBorder2{i}));
            end
            extendedMatrix1 = [matrix1(1,1) matrix1(1,:) matrix1(1,end);matrix1(:,1) matrix1 matrix1(:,end); matrix1(end,1) matrix1(end,:) matrix1(end,end)];
            extendedMatrix2 = [matrix2(1,1) matrix2(1,:) matrix2(1,end);matrix2(:,1) matrix2 matrix2(:,end); matrix2(end,1) matrix2(end,:) matrix2(end,end)];
        end
        
        function [matrixNext newTimeStep] = adaptiveTimeStep(cdag,matrixCurrent,changeTerm,oldTimeStep)
            matrixNext = matrixCurrent + cdag.dt*changeTerm;
            newTimeStep = oldTimeStep + 1;
            
            %take possibleDt as ratio of current state matrix and change in matrix,
            %align with default cdag.dt
            possibleDt = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
            possibleDt = possibleDt - mod(possibleDt,cdag.dt);
            if(possibleDt >= 1000*cdag.dt)
                possibleDt=1000*cdag.dt;
                matrixNext1 = matrixNext;
                
                [junk icdag.dx1] = min(abs(matrixCurrent./(matrixCurrent-matrixNext)));
                [junk icdag.dx2] = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
                %dimatrixp('HERE')
                %timematrixtep
                %icdag.dx1(icdag.dx2)
                %icdag.dx2
                %matrix(icdag.dx1(icdag.dx2),icdag.dx2)
                %matrixnext(icdag.dx1(icdag.dx2),icdag.dx2)
                %pomatrixmatrixibleDt
                %dimatrixp('HERE')
                %pomatrixmatrixibleDt
                
                %take step with current possibleDt, if relative change compared
                %with taking a step of size cdag.dt is greater than .01, halve
                %possibleDt
                matrixNext = matrix + possibleDt*changeTerm;
                %max(max(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnextmatrixtore)))
                while(max(max(abs((matrixNext-matrixNext1)./matrixNext1)))>.01)
                    possibleDt=possibleDt/2-mod(possibleDt/2,cdag.dt);
                    matrixNext = matrix + possibleDt*changeTerm;
                    %pomatrixmatrixibleDt
                end
                newTimeStep = oldTimeStep + uint64(possibleDt/cdag.dt);
                
                %pomatrixmatrixibleDt
                %max(max(matrix))
                %timematrixtep*cdag.dt
                %dimatrixp((max(max(matrix))-min(min(matrix)))/(max(max(matrix))))
                %dimatrixp((max(max(matrixnext))-min(min(matrixnext)))/(max(max(matrixnext))))
                %matrixubplot(ceil((length(matchematrixTimematrixToPrint)+1)/4),4,1+find(matchematrixTimematrixToPrint))
                %imagematrixc(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnext))
                %title(['t = ' num2matrixtr(uint64(timematrixtep*cdag.dt))])
                %aximatrix([1 matrixideLength 1 matrixideLength])
                %hold on
            end
        end
        
        function [cdag newOverflowCoord1 newOverflowCoord2 newOverflowXVal newOverflowLIKOVal newOverflowRIKOVal ...
                newOverflowLRKOVal newOverflowRRKOVal newOverflowIiVal newOverflowNiVal newOverflowColoniesVal] ...
                = overflowCell(varargin)
            
            cdag = varargin{1};
            
            %randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
            overflowCoord1 = varargin{2};
            overflowCoord2 = varargin{3};
            colony = varargin{4};
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
            
            newOverflowLIKOVal = cdag.LIKO(overflowCoord1,overflowCoord2);
            newOverflowRIKOVal = cdag.RIKO(overflowCoord1,overflowCoord2);
            newOverflowLRKOVal = cdag.LRKO(overflowCoord1,overflowCoord2);
            newOverflowRRKOVal = cdag.RRKO(overflowCoord1,overflowCoord2);
            newOverflowColoniesVal = cdag.colonies(overflowCoord1,overflowCoord2);
            if(first)
                newOverflowIiVal = cdag.Ii(overflowCoord1,overflowCoord2);
                newOverflowNiVal = cdag.Ni(overflowCoord1,overflowCoord2);
                newOverflowXVal = cdag.X(overflowCoord1,overflowCoord2)/2;
            else
                newOverflowIiVal = cdag.Ii(overflowCoord1,overflowCoord2);
                newOverflowNiVal = cdag.Ni(overflowCoord1,overflowCoord2);
                newOverflowXVal = cdag.X(overflowCoord1,overflowCoord2);
            end
            
            if(1)
                [XDist XNearestBorder] = bwdist(cdag.X==0 & cdag.onPlate,'chessboard');
                nearestBorderCoord1 = mod(XNearestBorder(overflowCoord1,overflowCoord2)-1,cdag.sideLength)+1;
                nearestBorderCoord2 = floor((double(XNearestBorder(overflowCoord1,overflowCoord2)) -1 )/ (cdag.sideLength*1.0))+1;
                randCoords=[];
                
                if(nearestBorderCoord1 < overflowCoord1)
                    randCoords(end+1,:) = [-1 0];
                end
                if(nearestBorderCoord1 > overflowCoord1)
                    randCoords(end+1,:) = [1 0];
                end
                if(nearestBorderCoord2 < overflowCoord2)
                    randCoords(end+1,:) = [0 -1];
                end
                if(nearestBorderCoord2 > overflowCoord2)
                    randCoords(end+1,:) = [0 1];
                end
                if(nearestBorderCoord1 < overflowCoord1 && nearestBorderCoord2 < overflowCoord2)
                    randCoords(end+1,:) = [-1 -1];
                end
                if(nearestBorderCoord1 < overflowCoord1 && nearestBorderCoord2 > overflowCoord2)
                    randCoords(end+1,:) = [-1 1];
                end
                if(nearestBorderCoord1 > overflowCoord1 && nearestBorderCoord2 < overflowCoord2)
                    randCoords(end+1,:) = [1 -1];
                end
                if(nearestBorderCoord1 > overflowCoord1 && nearestBorderCoord2 > overflowCoord2)
                    randCoords(end+1,:) = [1 1];
                end
                
                if(length(randCoords)==0)
                    newOverflowCoord1=overflowCoord1;
                    newOverflowCoord2=overflowCoord2;
                else
                    numRand = randi(size(randCoords,1),1);
                    newOverflowCoord1=overflowCoord1+randCoords(numRand,1);
                    if(newOverflowCoord1>cdag.sideLength)
                        newOverflowCoord1=1;
                    end
                    if(newOverflowCoord1<1)
                        newOverflowCoord1=cdag.sideLength;
                    end
                    newOverflowCoord2=overflowCoord2+randCoords(numRand,2);
                    if(newOverflowCoord2>cdag.sideLength)
                        newOverflowCoord2=1;
                    end
                    if(newOverflowCoord2<1)
                        newOverflowCoord2=cdag.sideLength;
                    end
                end
            end
            
            if(0)
                if(colony==0)
                    newOverflowCoord1 = overflowCoord1;
                    newOverflowCoord2 = overflowCoord2;
                elseif(overflowCoord1 == centerCoord1 && overflowCoord2 == centerCoord2)
                    randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
                    numRand = randi(size(randCoords,1),1);
                    newOverflowCoord1=overflowCoord1+randCoords(numRand,1);
                    if(newOverflowCoord1>cdag.sideLength)
                        newOverflowCoord1=1;
                    end
                    if(newOverflowCoord1<1)
                        newOverflowCoord1=cdag.sideLength;
                    end
                    newOverflowCoord2=overflowCoord2+randCoords(numRand,2);
                    if(newOverflowCoord2>cdag.sideLength)
                        newOverflowCoord2=1;
                    end
                    if(newOverflowCoord2<1)
                        newOverflowCoord2=cdag.sideLength;
                    end
                else
                    slope = (overflowCoord1 - centerCoord1) / (overflowCoord2 - centerCoord2);
                    if(isnan(slope))
                        if(overflowCoord1<centerCoord1)
                            newOverflowCoord1=overflowCoord1-1;newOverflowCoord2=overflowCoord2;
                        end
                        if(overflowCoord1>centerCoord1)
                            newOverflowCoord1=overflowCoord1+1;newOverflowCoord2=overflowCoord2;
                        end
                    else
                        if(slope<=1/2 && slope>-1/2)
                            if(overflowCoord2<centerCoord2)
                                newOverflowCoord1=overflowCoord1;newOverflowCoord2=overflowCoord2-1;
                            end
                            if(overflowCoord2>centerCoord2)
                                newOverflowCoord1=overflowCoord1;newOverflowCoord2=overflowCoord2+1;
                            end
                        elseif(slope<=-1/2 && slope>-2)
                            if(overflowCoord2<centerCoord2)
                                newOverflowCoord1=overflowCoord1-1;newOverflowCoord2=overflowCoord2-1;
                            end
                            if(overflowCoord2>centerCoord2)
                                newOverflowCoord1=overflowCoord1+1;newOverflowCoord2=overflowCoord2+1;
                            end
                        elseif(slope<=-2 || slope>2)
                            if(overflowCoord1<centerCoord1)
                                newOverflowCoord1=overflowCoord1-1;newOverflowCoord2=overflowCoord2;
                            end
                            if(overflowCoord1>centerCoord1)
                                newOverflowCoord1=overflowCoord1+1;newOverflowCoord2=overflowCoord2;
                            end
                        elseif(slope<=2 && slope > 1/2)
                            if(overflowCoord2<centerCoord2)
                                newOverflowCoord1=overflowCoord1+1;newOverflowCoord2=overflowCoord2-1;
                            end
                            if(overflowCoord2>centerCoord2)
                                newOverflowCoord1=overflowCoord1-1;newOverflowCoord2=overflowCoord2+1;
                            end
                        end
                    end
                    if(newOverflowCoord1>cdag.sideLength)
                        newOverflowCoord1=1;
                    end
                    if(newOverflowCoord1<1)
                        newOverflowCoord1=cdag.sideLength;
                    end
                    if(newOverflowCoord2>cdag.sideLength)
                        newOverflowCoord2=1;
                    end
                    if(newOverflowCoord2<1)
                        newOverflowCoord2=cdag.sideLength;
                    end
                end
            end
            
            if(first)
                cdag.Ii(overflowCoord1,overflowCoord2) = cdag.Ii(overflowCoord1,overflowCoord2);
                cdag.Ni(overflowCoord1,overflowCoord2) = cdag.Ni(overflowCoord1,overflowCoord2);
                cdag.X(overflowCoord1,overflowCoord2) = cdag.X(overflowCoord1,overflowCoord2)/2;
            else
                cdag.LIKO(overflowCoord1,overflowCoord2) = overflowLIKOVal;
                cdag.RIKO(overflowCoord1,overflowCoord2) = overflowRIKOVal;
                cdag.LRKO(overflowCoord1,overflowCoord2) = overflowLRKOVal;
                cdag.RRKO(overflowCoord1,overflowCoord2) = overflowRRKOVal;
                cdag.colonies(overflowCoord1,overflowCoord2) = overflowColoniesVal;
                cdag.Ii(overflowCoord1,overflowCoord2) = overflowIiVal;
                cdag.Ni(overflowCoord1,overflowCoord2) = overflowNiVal;
                cdag.X(overflowCoord1,overflowCoord2) = overflowXVal;
            end
        end
        
        function main(cdag)
            tic
            while(cdag.timeStep <= cdag.numTimeSteps)
                cdag.runOneTimeStep();
            end
            toc
            cdag.makeGraphics();
        end
        
        function cdag = runOneTimeStep(cdag)
            %disp(cdag.X)
            disp(num2str(cdag.timeStep*cdag.dt))
            %disp(num2str(max(max(cdag.Ii))))
            
            % cdag.L cdag.dynamics
            diffusionTermL = cdag.makeDiffusionMatrix(cdag.L,cdag.DL,cdag.extendedL1,cdag.extendedL2);
            %diffusionTermL = cdag.DL*filter2(diffKernel,extendedL,'valid')/(cdag.dx)^2;
            productionTermL = cdag.kprodL*(cdag.X.*(~cdag.LIKO));
            if(cdag.withAutoInduction)
                productionTermL = productionTermL + cdag.kprodL*(cdag.X.*(~cdag.LIKO).*(~cdag.LRKO)).*cdag.L./(cdag.kautoL+cdag.L);
            end
            decayTermL = cdag.kdecayL*cdag.L;
            changeTermL = diffusionTermL + productionTermL - decayTermL;
            
            % cdag.R cdag.dynamics
            diffusionTermR = cdag.makeDiffusionMatrix(cdag.R,cdag.DR,cdag.extendedR1,cdag.extendedR2);
            %diffusionTermR = DS*filter2(diffKernel,extendedR,'valid')/(cdag.dx)^2;
            productionTermR = cdag.kprodR*(cdag.X.*(~cdag.RIKO));
            if(cdag.withAutoInduction)
                productionTermR = productionTermR + (cdag.kprodR*(cdag.X.*(~cdag.RIKO).*(~cdag.RRKO).*(~cdag.LRKO).*cdag.R./(cdag.kautoR+cdag.R)).* (cdag.L./(cdag.kcrossL + cdag.L)));
            end
            decayTermR = cdag.kdecayR*cdag.R;
            changeTermR = diffusionTermR + productionTermR - decayTermR;
            
            diffusionTermC = cdag.makeDiffusionMatrix(cdag.C,cdag.DC,cdag.extendedC1,cdag.extendedC2);
            diffusionTermN = cdag.makeDiffusionMatrix(cdag.N,cdag.DN,cdag.extendedN1,cdag.extendedN2);
            diffusionTermI = cdag.makeDiffusionMatrix(cdag.I,cdag.DI,cdag.extendedI1,cdag.extendedI2);
            changeTermC = diffusionTermC;
            changeTermN = diffusionTermN;
            changeTermI = diffusionTermI;
            cdag.C = max(0,cdag.C + cdag.dt*changeTermC);
            cdag.N = max(0,cdag.N + cdag.dt*changeTermN);
            cdag.I = max(0,cdag.I + cdag.dt*changeTermI);
            cdag.mu = cdag.muMax*(cdag.C>0 & cdag.N>0).*cdag.Ii + cdag.muMaxPrime*(cdag.N<=0 & cdag.C>0).*cdag.Ni;
            consumptionTermC = -1/cdag.YC*cdag.mu.*cdag.X.*(cdag.C>0);
            consumptionTermN = -1/cdag.YN*cdag.mu.*cdag.X.*(cdag.N>0);
            consumptionTermI = -1/cdag.YI*cdag.mu.*cdag.X.*(cdag.I>0);
            consumptionTermNi = -(1/cdag.YNi+cdag.Ni).*cdag.mu.*(cdag.N==0 & cdag.Ni>0);
            consumptionTermIi = -(1/cdag.YIi+cdag.Ii).*cdag.mu.*(cdag.I==0 & cdag.Ii>0);
            growthTermX = cdag.mu.*cdag.X.*(cdag.C>0) + (cdag.mu-cdag.kdecayX).*cdag.X.*(cdag.C==0);
            %disp(num2str(max(max(growthTermX))))
            changeTermC = consumptionTermC;
            changeTermN = consumptionTermN;
            changeTermI = consumptionTermI;
            changeTermNi = consumptionTermNi;
            changeTermIi = consumptionTermIi;
            changeTermX = growthTermX;
            cdag.C = max(0,cdag.C + cdag.dt*changeTermC);
            cdag.N = max(0,cdag.N + cdag.dt*changeTermN);
            cdag.I = max(0,cdag.I + cdag.dt*changeTermI);
            cdag.Ni = max(0,cdag.Ni + cdag.dt*changeTermNi);
            cdag.Ii = max(0,cdag.Ii + cdag.dt*changeTermIi);
            cdag.X = max(0,cdag.X + cdag.dt*changeTermX);
            if(min(min(consumptionTermC))~=0)
                cdag.C(cdag.C<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermC(consumptionTermC~=0))))=0;
            end
            if(min(min(consumptionTermN))~=0)
                %-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermN(consumptionTermN~=0)))
                cdag.N(cdag.N<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermN(consumptionTermN~=0))))=0;
            end
            if(min(min(consumptionTermI))~=0)
                cdag.I(cdag.I<-1/(cdag.sideLength*cdag.sideLength)*max(max(cdag.dt*consumptionTermI(consumptionTermI~=0))))=0;
            end
            
            [overflowCoord1 overflowCoord2] = find(cdag.X>cdag.XThresh,1,'first');
            %randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
            while(~isempty(overflowCoord1))
                [cdag overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
                    overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal] = ...
                    cdag.overflowCell(overflowCoord1,overflowCoord2, ...
                    cdag.colonies(overflowCoord1,overflowCoord2));
                while( overflowXVal > 0 )
                    [cdag overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
                        overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal] = ...
                        cdag.overflowCell(overflowCoord1,overflowCoord2, ...
                        cdag.colonies(overflowCoord1,overflowCoord2), ...
                        overflowXVal,overflowLIKOVal,overflowRIKOVal,overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal,overflowColoniesVal);
                end
                [overflowCoord1 overflowCoord2] = find(cdag.X>cdag.XThresh,1,'first');
            end
            
            %implicit equation
            %forms without decay and with autoinduction?
            %Rnext = (cdag.R+cdag.dt*...
            %(kprod*cdag.X+DS*([cdag.R(:,2:end) cdag.R(:,1)] + [cdag.R(:,end) cdag.R(:,1:end-1)] - ...
            %4*cdag.R + [cdag.R(2:end, :); cdag.R(1,:)] + [cdag.R(end,:); cdag.R(1:end-1, :)])...
            %/(cdag.dx)^2)/(1+cdag.dt*kdecay));
            
            if(mod(cdag.timeStep,uint64(100/cdag.dt))==0)
                %         disp(cdag.dt*cdag.timeStep)
                %         disp(max(max(cdag.R)))
                %         disp((max(max(cdag.R))-min(min(cdag.R)))/(max(max(cdag.R))))
                %         disp(sum(sum(decayTermR)))
                %         disp(sum(sum(productionTermR)))
                %         cdag.Ii(cdag.quarterPt,cdag.quarterPt)
                %         cdag.Ii(cdag.quarterPt,cdag.quarterPt*3)
                %         cdag.Ii(cdag.quarterPt*3,cdag.quarterPt*3)
                %         cdag.Ii(cdag.quarterPt*3,cdag.quarterPt)
            end
            
            cdag.xCum(cdag.timeStep) = sum(sum(cdag.X));
            cdag.cCum(cdag.timeStep) = sum(sum(cdag.C));
            cdag.nCum(cdag.timeStep) = sum(sum(cdag.N));
            cdag.iCum(cdag.timeStep) = sum(sum(cdag.I));
            cdag.rCum(cdag.timeStep) = sum(sum(cdag.R));
            cdag.lCum(cdag.timeStep) = sum(sum(cdag.L));
            
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
            
            if(cdag.withAdaptive)
                [Lnext newTimeStepL] = adaptiveTimeStep(cdag.L,cdag.dt,changeTermL,cdag.timeStep);
                [Rnext newTimeStepR] = adaptiveTimeStep(cdag.R,cdag.dt,changeTermR,cdag.timeStep);
                Lnext = cdag.L + min(newTimeStepL, newTimeStepR) * (changeTermL);
                Rnext = cdag.R + min(newTimeStepL, newTimeStepR) * (changeTermR);
                %cdag.timeStep = min(new
            else
                Lnext = cdag.L + cdag.dt*changeTermL;
                Rnext = cdag.R + cdag.dt*changeTermR;
                GFPnext = cdag.R.*cdag.X.*cdag.RIKO;
                cdag.timeStep = cdag.timeStep + 1;
            end
            
            cdag.Lstore(cdag.timeStep) = cdag.L(cdag.midPt,cdag.midPt);
            cdag.L = Lnext;
            %extendedL = [cdag.L(end,end) cdag.L(end,:) cdag.L(end,1);cdag.L(:,end) cdag.L cdag.L(:,1); cdag.L(1,end) cdag.L(1,:) cdag.L(1,1)];
            
            cdag.Rstore(cdag.timeStep) = cdag.R(cdag.midPt,cdag.midPt);
            cdag.R = Rnext;
            %extendedR = [cdag.R(end,end) cdag.R(end,:) cdag.R(end,1);cdag.R(:,end) cdag.R cdag.R(:,1); cdag.R(1,end) cdag.R(1,:) cdag.R(1,1)];
            
            cdag.GFP = GFPnext;
            
            %     extendedL = [cdag.L(1,1) cdag.L(1,:) cdag.L(1,end);cdag.L(:,1) cdag.L cdag.L(:,end); cdag.L(end,1) cdag.L(end,:) cdag.L(end,end)];
            %     extendedR = [cdag.R(1,1) cdag.R(1,:) cdag.R(1,end);cdag.R(:,1) cdag.R cdag.R(:,end); cdag.R(end,1) cdag.R(end,:) cdag.R(end,end)];
            %     extendedC = [cdag.C(1,1) cdag.C(1,:) cdag.C(1,end);cdag.C(:,1) cdag.C cdag.C(:,end); cdag.C(end,1) cdag.C(end,:) cdag.C(end,end)];
            %     extendedN = [cdag.N(1,1) cdag.N(1,:) cdag.N(1,end);cdag.N(:,1) cdag.N cdag.N(:,end); cdag.N(end,1) cdag.N(end,:) cdag.N(end,end)];
            %     extendedI = [cdag.I(1,1) cdag.I(1,:) cdag.I(1,end);cdag.I(:,1) cdag.I cdag.I(:,end); cdag.I(end,1) cdag.I(end,:) cdag.I(end,end)];
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
                        saveMap=matrixToPlotScale;
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
            
            %if(0)
            system('rm frames/*2*.jpg')
            MINVAL1=min(min(min(cdag.GFPAll2)));
            MAXVAL1=max(max(max(cdag.GFPAll2)));
            MINVAL2=min(min(min(cdag.XAll2)));
            MAXVAL2=max(max(max(cdag.XAll2)));
            MINVAL3=min(min(min(cdag.RAll2)));
            MAXVAL3=max(max(max(cdag.RAll2)));
            for i=1:length(cdag.frameTimeSteps)
                cdag.GFP = cdag.GFPAll2(:,:,i);
                cdag.X = cdag.XAll2(:,:,i);
                cdag.RIKO = cdag.RIKOAll2(:,:,i);
                cdag.R = cdag.RAll2(:,:,i);
                %cdag.GFP = cdag.GFP + cdag.RIKO*.1*MAXVAL1;
                matrixToPlot(:, :, 1) = cdag.X.*(~cdag.RIKO);
                matrixToPlot(:, :, 2) = cdag.GFP;
                matrixToPlot(:, :, 3) = cdag.R.*(cdag.X==0);
                matrixToPlotScale = uint16(zeros(size(matrixToPlot,1),size(matrixToPlot,2),size(matrixToPlot,3)));
                matrixToPlotScale(:, :, 3) = uint16( (matrixToPlot(:, :, 3) - MINVAL3) / (MAXVAL3 - MINVAL3) * (2^16-1) );
                matrixToPlotScale(:, :, 2) = uint16( (matrixToPlot(:, :, 2) - MINVAL1) / (MAXVAL1 - MINVAL1) * (2^16-1) );
                matrixToPlotScale(:, :, 2) = min( uint16(matrixToPlotScale(:, :, 2) + uint16(cdag.RIKO*2000)), uint16(2^16-1));
                matrixToPlotScale(:, :, 1) = uint16( (matrixToPlot(:, :, 1) - MINVAL2) / (MAXVAL2 - MINVAL2) * (2^16-1) );
                figure(15)
                imshow(matrixToPlotScale(cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,:))
                drawnow
                print( 15, '-djpeg', sprintf('frames/frmGFP4%d.jpg', i))
            end
            
            disp('making ffmpeg movie')
            system('"cdag.C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmGFP4%d.jpg frames/GFP4.mov')
            disp('done with ffmpeg movie')
            
            MINVAL1=min(min(min(cdag.CAll2)));
            MAXVAL1=max(max(max(cdag.CAll2)));
            MINVAL2=min(min(min(cdag.XAll2)));
            MAXVAL2=max(max(max(cdag.XAll2)));
            for i=1:length(cdag.frameTimeSteps)
                cdag.C = cdag.CAll2(:,:,i);
                cdag.X = cdag.XAll2(:,:,i);
                matrixToPlot(:, :, 1) = cdag.X;
                matrixToPlot(:, :, 2) = uint16(0);
                matrixToPlot(:, :, 3) = cdag.C;
                matrixToPlotScale = uint16(zeros(size(matrixToPlot,1),size(matrixToPlot,2),size(matrixToPlot,3)));
                matrixToPlotScale(:,:,3) = uint16( (matrixToPlot(:,:,3) - MINVAL1) / (MAXVAL1 - MINVAL1) * (2^16-1) );
                matrixToPlotScale(:,:,1) = uint16( (matrixToPlot(:,:,1) - MINVAL2) / (MAXVAL2 - MINVAL2) * (2^16-1) );
                figure(15)
                imshow(matrixToPlotScale(cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,cdag.sideOffset+1:cdag.sideLength-cdag.sideOffset,:))
                drawnow
                print( 15, '-djpeg', sprintf('frames/frmX4%d.jpg', i))
            end
            
            disp('making ffmpeg movie')
            system('"cdag.C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmX4%d.jpg frames/X4.mov')
            disp('done with ffmpeg movie')
            %end
        end
    end
end