classdef coloniesDiffusionAndGrowth
    properties (SetAccess = protected)
        withAutoInduction = 1;
        withAdaptive = 0;
        
        sideLength = 200;
        midPt = ceil(sideLength/2.);
        quarterPt = floor(midPoint/2);
        
        DR = .030; % Diffusivity mm^2/min
        DL = .030;
        DC = .0556; %bionumbers
        DN = .045; %http://www.sciencedirect.com/science/article/pii/S0022024800005285
        DI = .045; %assumed to be close to DN
        
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
        dt = dx^2/max([DR,DL,DC,DN,DI])*.2495; %CFL constant
        tend = 2400; % End time (s)
        % timesToPrint = [1 25 50 75 100 200 300 500 1000];
        timesToPrint = [1 250 500 750 1000 1250 1500 1750 2000];
        %timesToPrint = dt*[1 2 3 4 5 6 7 8 9 10];
        cmap = jet(length(timesToPrint));
        numTimeSteps = length(0:dt:tend);
        frameTimeSteps = 1:ceil(numTimeSteps/100):numTimeSteps;
        
        X = zeros(sideLength, sideLength);
        XThresh = 10000*1/(.15*1.1*10^9);%.0001
        LIKO = zeros(sideLength, sideLength);
        RIKO = zeros(sideLength, sideLength);
        LRKO = zeros(sideLength, sideLength);
        RRKO = zeros(sideLength, sideLength);
        colonies = zeros(sideLength, sideLength);
        R = zeros(sideLength, sideLength);
        L = zeros(sideLength, sideLength);
        initialR = 0; %.0001;
        initialL = 0;
        R(midPoint,midPoint) = initialR;
        L(midPoint,midPoint) = initialR;
        Rstore = zeros(1,numTimeSteps);
        Lstore = zeros(1,numTimeSteps);
        GFP = R.*X.*RIKO;
        C = .5*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
        N = .0625*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
        I = .00006981*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
        Ni = ones(sideLength, sideLength).*(X~=0);
        Ii = ones(sideLength, sideLength).*(X~=0);
        mu = zeros(sideLength, sideLength);
        muMax = 0.3341/60; % hourly rate scaled to minutes
        muMaxPrime = .0616/60;
        
        sideOffset=0;
        X(ceil(sideOffset+(sideLength-2*sideOffset)*100/350),ceil(sideOffset+(sideLength-2*sideOffset)*50/350))=100/(.15*1.1*10^11);
        colonies(ceil(sideOffset+(sideLength-2*sideOffset)*100/350),ceil(sideOffset+(sideLength-2*sideOffset)*50/350))=1;
        colonyCenterCoords = [ceil(sideOffset+(sideLength-2*sideOffset)*100/350) ceil(sideOffset+(sideLength-2*sideOffset)*50/350)];
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
            X(ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350),ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350))=100/(.15*1.1*10^11);
            RIKO(ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350),ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350))=1;
            colonies(ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350),ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350))=i+1;
            colonyCenterCoords(end+1,:)=[ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350) ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350)];
        end
        
        RAll = zeros(sideLength, sideLength, length(timesToPrint));
        LAll = zeros(sideLength, sideLength, length(timesToPrint));
        XAll = zeros(sideLength, sideLength, length(timesToPrint));
        GFPAll = zeros(sideLength, sideLength, length(timesToPrint));
        GFPAll2 = zeros(sideLength, sideLength, length(frameTimeSteps));
        XAll2 = zeros(sideLength, sideLength, length(frameTimeSteps));
        CAll2 = zeros(sideLength, sideLength, length(frameTimeSteps));
        RIKOAll2 = zeros(sideLength, sideLength, length(frameTimeSteps));
        RAll2 = zeros(sideLength, sideLength, length(frameTimeSteps));
        
        [coords1 coords2] = meshgrid(1:sideLength,1:sideLength);
        onPlate = (sqrt((coords1-midPoint).^2 + (coords2-midPoint).^2) <= (midPoint-1));
        horzPlateEdge = ~onPlate & ([zeros(sideLength,1) onPlate(:,1:end-1)] | [onPlate(:,2:end) zeros(sideLength,1)]);
        vertPlateEdge = ~onPlate & ([zeros(1,sideLength); onPlate(1:end-1,:)] | [onPlate(2:end,:); zeros(1,sideLength)]);
        pointsBorderToInside1 = containers.Map('KeyType','double','ValueType','double');
        pointsBorderToInside2 = containers.Map('KeyType','double','ValueType','double');
        for i=1:sideLength
            if(abs(diff(find(horzPlateEdge(i,:))))==2)
                jCoords=find(horzPlateEdge(i,:));
                pointsBorderToInside1((jCoords(1)-1)*sideLength+i)=jCoords(1)*sideLength+i;
                pointsBorderToInside1((jCoords(2)-1)*sideLength+i)=jCoords(1)*sideLength+i;
                %else
                horzPlateEdge(i,:) = zeros(1,sideLength);
            end
        end
        for i=1:sideLength
            if(abs(diff(find(vertPlateEdge(:,i))))==2)
                jCoords=find(vertPlateEdge(:,i));
                pointsBorderToInside2((i-1)*sideLength+jCoords(1))=(i-1)*sideLength+jCoords(1)+1;
                pointsBorderToInside2((i-1)*sideLength+jCoords(2))=(i-1)*sideLength+jCoords(1)+1;
                %else
                vertPlateEdge(:,i) = zeros(sideLength,1);
            end
        end
        
        diffKernel1 = [0 0 0; 1 -2 1; 0 0 0];
        diffKernel2 = [0 1 0; 0 -2 0; 0 1 0];
        
        X(~onPlate) = 0;
        R(~onPlate) = 0;
        L(~onPlate) = 0;
        LIKO(~onPlate) = 0;
        RIKO(~onPlate) = 0;
        LRKO(~onPlate) = 0;
        RRKO(~onPlate) = 0;
        C(~onPlate) = 0;
        N(~onPlate) = 0;
        I(~onPlate) = 0;
        Ni(~onPlate) = 0;
        Ii(~onPlate) = 0;
        [extendedL1 extendedL2] = makeExtendedMatrix(L);
        [extendedR1 extendedR2] = makeExtendedMatrix(R);
        [extendedC1 extendedC2] = makeExtendedMatrix(C);
        [extendedN1 extendedN2] = makeExtendedMatrix(N);
        [extendedI1 extendedI2] = makeExtendedMatrix(I);
    end
end