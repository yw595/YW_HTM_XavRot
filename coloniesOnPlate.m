% This is a pde model for quorum sensing factor diffusion, production, and
% decay
% 
%
% Authors: Hilary Monaco and Yiping Wang
% 
% Date Updated: July 25, 2014

clear all; close all; clc

%%% Options %%%

withAutoInduction = 1;
withAdaptive = 0;

%%% Constants %%%

sideLength = 20;
midPoint = ceil(sideLength/2.);
quarterPt = floor(midPoint/2);

DR = .010; % Diffusivity mm^2/min
DL = .010; 
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
kdecayR = .0001; %Decay rate molecules/second
kdecayL = .0001; %Decay rate molecules/second
kdecayX = .0062/60; %hourly rate scaled to minutes
kautoR = 0.0001;
kautoL = 0.0001;
kcrossL = 0.0005;

muMax = 0.3341/60; % hourly rate scaled to minutes
muMaxPrime = .0616/60;

dx = .1; %units of mm
dy = .1; %units of mm
dt = dx^2/max([DR,DL,DC,DN,DI])*.2495; %CFL constant 
tend = 1000; % End time (s)

XThresh = 2;

timesToPrint = [1 250 500 750 1000];
cmap = jet(length(timesToPrint));
numTimeSteps = length(0:dt:tend);
frameTimeSteps = 1:ceil(numTimeSteps/100):numTimeSteps;
Rstore = zeros(1,numTimeSteps);
Lstore = zeros(1,numTimeSteps);

%%% Initial Conditions %%%

X = zeros(sideLength, sideLength); 
LIKO = zeros(sideLength, sideLength); 
RIKO = zeros(sideLength, sideLength);
LRKO = zeros(sideLength, sideLength); 
RRKO = zeros(sideLength, sideLength); 
%XDKO = zeros(sideLength, sideLength); 

% %%% Multiple Colonies
X(midPoint,midPoint) = 1;
X(quarterPt, quarterPt) = 1;
X(quarterPt, quarterPt*3) = 1;
X(quarterPt*3, quarterPt*3) = 1;
X(quarterPt*3, quarterPt) = 1;
LIKO(quarterPt, quarterPt) = 1;
RIKO(quarterPt, quarterPt*3) = 1;
LRKO(quarterPt*3, quarterPt*3) = 1;
RRKO(quarterPt*3, quarterPt) = 1;

sideOffset=0;
% X(ceil(sideOffset+(sideLength-2*sideOffset)*100/350),ceil(sideOffset+(sideLength-2*sideOffset)*50/350))=1;
% RIKOcoords=[160 25;
%     210 25;
%     275 65;
%     110 75;
%     80 100;
%     310 100;
%     40 120;
%     125 120;
%     175 120;
%     35 135;
%     280 135;
%     145 150;
%     310 150;
%     20 160;
%     260 170;
%     90 190;
%     220 200;
%     325 210;
%     55 220;
%     105 225;
%     125 250;
%     155 275;
%     120 290;
%     220 295;
%     160 300];
% for i=1:length(RIKOcoords)
%     X(ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350),ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350))=1;
%     RIKO(ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,1)/350),ceil(sideOffset+(sideLength-2*sideOffset)*RIKOcoords(i,2)/350))=1;
% end

R = zeros(sideLength, sideLength);
L = zeros(sideLength, sideLength);
GFP = R.*X.*RIKO;
C = .5*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
N = .0625*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
%P = .01*ones(sideLength, sideLength);
I = .00006981*.0012*1/(sideLength*sideLength)*ones(sideLength, sideLength);
Ni = ones(sideLength, sideLength).*(X~=0);
Ii = ones(sideLength, sideLength).*(X~=0);
mu = zeros(sideLength, sideLength);
RAll = zeros(sideLength, sideLength, length(timesToPrint));
LAll = zeros(sideLength, sideLength, length(timesToPrint));
XAll = zeros(sideLength, sideLength, length(timesToPrint));
GFPAll = zeros(sideLength, sideLength, length(timesToPrint));
RIKOAll = zeros(sideLength, sideLength, length(frameTimeSteps));
initialR = 0; %.0001;
initialL = 0;
R(midPoint,midPoint) = initialR;
L(midPoint,midPoint) = initialR;

%%% Solving PDE Model %%%

%diffKernel = [0 1 0; 1 -4 1; 0 1 0];
%extendedL = [L(end,end) L(end,:) L(end,1);L(:,end) L L(:,1); L(1,end) L(1,:) L(1,1)];
%extendedR = [R(end,end) R(end,:) R(end,1);R(:,end) R R(:,1); R(1,end) R(1,:) R(1,1)];
% This is the the explicit and implicit methods
skipsteps=0;
tic
timeStep=1;
while(timeStep <= numTimeSteps)
    mu = muMax*(C>0 & N>0).*Ii + muMaxPrime*(N<=0 & C>0).*Ni;
    % timeStep*dt
    
    % L dynamics
    diffusionTermL = DL*([L(:,2:end) L(:,1)] + [L(:,end) L(:,1:end-1)] - ...
        4*L + [L(2:end, :); L(1,:)] + [L(end,:); L(1:end-1, :)])/(dx)^2;
    %diffusionTermL = DL*filter2(diffKernel,extendedL,'valid')/(dx)^2;
    productionTermL = kprodL*(X.*(~LIKO));
    if(withAutoInduction)
        productionTermL = productionTermL + kprodL*(X.*(~LIKO).*(~LRKO)).*L./(kautoL+L);
    end
    decayTermL = kdecayL*L;
    changeTermL = diffusionTermL + productionTermL - decayTermL;
    
    % R dynamics
    diffusionTermR = DR*([R(:,2:end) R(:,1)] + [R(:,end) R(:,1:end-1)] - ...
        4*R + [R(2:end, :); R(1,:)] + [R(end,:); R(1:end-1, :)])/(dx)^2;
    %diffusionTermR = DS*filter2(diffKernel,extendedR,'valid')/(dx)^2;
    productionTermR = kprodR*(X.*(~RIKO));
    if(withAutoInduction)
        productionTermR = productionTermR + (kprodR*(X.*(~RIKO).*(~RRKO).*(~LRKO).*R./(kautoR+R)).* (L./(kcrossL + L)));
    end
    decayTermR = kdecayR*R;
    changeTermR = diffusionTermR + productionTermR - decayTermR;
    
    diffusionTermC = DC*([C(:,2:end) C(:,1)] + [C(:,end) C(:,1:end-1)] - ...
        4*C + [C(2:end, :); C(1,:)] + [C(end,:); C(1:end-1, :)])/(dx)^2;
    diffusionTermN = DN*([N(:,2:end) N(:,1)] + [N(:,end) N(:,1:end-1)] - ...
        4*N + [N(2:end, :); N(1,:)] + [N(end,:); N(1:end-1, :)])/(dx)^2;
    diffusionTermI = DI*([I(:,2:end) I(:,1)] + [I(:,end) I(:,1:end-1)] - ...
        4*I + [I(2:end, :); I(1,:)] + [I(end,:); I(1:end-1, :)])/(dx)^2;
    consumptionTermC = -1/YC*mu.*X.*(C>0);
    consumptionTermN = -1/YN*mu.*X.*(N>0);
    consumptionTermI = -1/YI*mu.*X.*(I>0);
    consumptionTermNi = -(1/YNi + Ni).*mu.*(N==0 & Ni>0);
    consumptionTermIi = -(1/YIi + Ii).*mu.*(I==0 & Ii>0);
    growthTermX = mu.*X.*(C>0) + (mu-kdecayX).*X.*(C==0);
    changeTermC = diffusionTermC + consumptionTermC;
    changeTermN = diffusionTermN + consumptionTermN;
    changeTermI = diffusionTermI + consumptionTermI;
    changeTermNi = consumptionTermNi;
    changeTermIi = consumptionTermIi;
    changeTermX = growthTermX;
    C = max(0,C + dt*changeTermC);
    N = max(0,N + dt*changeTermN);
    I = max(0,I + dt*changeTermI);
    Ni = max(0,Ni + dt*changeTermNi);
    Ii = max(0,Ii + dt*changeTermIi);
    X = max(0,X + dt*changeTermX);
    
    [overflowCoord1 overflowCoord2] = find(X>XThresh,1,'first');
    %randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
    while(~isempty(overflowCoord1))
        [LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
        overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal] = ...
        overflowCell(LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,overflowCoord1,overflowCoord2);
        while( overflowXVal > 0 )
            [LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
            overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal] = ...
            overflowCell(LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
            overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal);
        end
        [overflowCoord1 overflowCoord2] = find(X>XThresh,1,'first');
    end
    
    %implicit equation
    %forms without decay and with autoinduction?
    %Rnext = (R+dt*...
    %(kprod*X+DS*([R(:,2:end) R(:,1)] + [R(:,end) R(:,1:end-1)] - ...
    %4*R + [R(2:end, :); R(1,:)] + [R(end,:); R(1:end-1, :)])...
    %/(dx)^2)/(1+dt*kdecay));
    
    if(mod(timeStep,uint64(100/dt))==0)
        disp(dt*timeStep)
        disp(max(max(R)))
        disp((max(max(R))-min(min(R)))/(max(max(R))))
        disp(sum(sum(decayTermR)))
        disp(sum(sum(productionTermR)))
        Ii(quarterPt,quarterPt)
        Ii(quarterPt,quarterPt*3)
        Ii(quarterPt*3,quarterPt*3)
        Ii(quarterPt*3,quarterPt)
    end
    
    xCum(timeStep) = sum(sum(X));
    cCum(timeStep) = sum(sum(C));
    nCum(timeStep) = sum(sum(N));
    iCum(timeStep) = sum(sum(I));
    rCum(timeStep) = sum(sum(R));
    lCum(timeStep) = sum(sum(L));
    
    matchesTimesToPrint = (uint64(timesToPrint/dt) == timeStep);
    if sum(matchesTimesToPrint)~=0
        RAll(:,:,matchesTimesToPrint) = R;
        LAll(:,:,matchesTimesToPrint) = L;
        XAll(:,:,matchesTimesToPrint) = X;
        GFPAll(:,:,matchesTimesToPrint) = GFP;
        RIKOAll(:,:,frameTimeSteps == timeStep) = RIKO;
    end
    
    if(withAdaptive)
        [Lnext newTimeStepL] = adaptiveTimeStep(L,dt,changeTermL,timeStep);
        [Rnext newTimeStepR] = adaptiveTimeStep(R,dt,changeTermR,timeStep);
        Lnext = L + min(newTimeStepL, newTimeStepR) * (changeTermL);
        Rnext = R + min(newTimeStepL, newTimeStepR) * (changeTermR);
        %timeStep = min(new
    else
        Lnext = L + dt*changeTermL;
        Rnext = R + dt*changeTermR;
        GFPnext = R.*X.*RIKO;
        timeStep = timeStep + 1;
    end
    
    Lstore(timeStep) = L(midPoint,midPoint);
    L = Lnext;
    %extendedL = [L(end,end) L(end,:) L(end,1);L(:,end) L L(:,1); L(1,end) L(1,:) L(1,1)];
    
    Rstore(timeStep) = R(midPoint,midPoint);
    R = Rnext;
    %extendedR = [R(end,end) R(end,:) R(end,1);R(:,end) R R(:,1); R(1,end) R(1,:) R(1,1)];
    
    GFP = GFPnext;
    
end
toc

makeMovie=0;
titles={'L','R','L and R','X','GFP'};
for i=1:length(titles)
    if strcmp(titles{i},'L')
        arrayAll = LAll;
    elseif strcmp(titles{i},'R')
        arrayAll = RAll;
    elseif strcmp(titles{i},'L and R')
        arrayAll = [RAll LAll];
    elseif strcmp(titles{i},'X')
        arrayAll = [XAll CAll];
    elseif strcmp(titles{i},'GFP')
        arrayAll = [GFPAll XAll RIKOAll];
    end
    makeGraphics(arrayAll,timesToPrint,i,titles{i},cmap,makeMovie);
end

figure(6)
surf(C)
figure(7)
surf(N)
figure(8)
surf(I)
figure(9)
surf(Ni)
figure(11)
surf(Ii)

% figure(11)
% hold on
% subplot(3,4,1)
% plot(1:numTimeSteps, xCum, 'k')
% subplot(3,4,2)
% plot(1:numTimeSteps, cCum, 'm')
% subplot(3,4,3)
% plot(1:numTimeSteps, nCum, 'c')
% subplot(3,4,4)
% plot(1:numTimeSteps, iCum, 'y')
% subplot(3,4,5)
% plot(1:numTimeSteps, rCum, 'g')
% subplot(3,4,6)
% imshow(C)
% subplot(3,4,7)
% imshow(N)
% subplot(3,4,8)
% imshow(I)
% subplot(3,4,9)
% plot(1:numTimeSteps, lCum, 'b')

% This is the Matrix method without autoinduction
if(0)
    
    R = reshape(R,sideLength^2,1)+.01;
    X = reshape(X,sideLength^2,1);
    % reshape(R, M, N) reshapes R -> MXN matrix
    
    %multiplicative induction
    %kprod=kdecay;
    
    b = kprodR*X; % Density Dependent production
    %b=zeros(sideLength^2,1);
    
    % Create operator matrix
    bigEye = eye(sideLength^2); % creates basline identity matrix
    
    A = DS/(dx^2)*(zeros(sideLength^2,sideLength^2)-4*bigEye+...
        circshift(bigEye,[0 1])+circshift(bigEye,[0 -1])+...
        circshift(bigEye,[0 sideLength])+circshift(bigEye,[0 -sideLength]));
    % circshift(A, [y X]) shifts the values of the array by y vertically
    % and X horizontally
    
    A = A-kdecayR*bigEye;
    
    %multiplicative induction
    %A=A-kdecay*bigEye+diag(kprod*X);
    
    A = sparse(A);
    
    % Rolve for R (Concentration of R everywhere)
    R = A\(-b); % R = Del^2\(-production) at steady state
    %max(abs(A*R+b))
    
    R = reshape(R,sideLength,sideLength);
    
    %figure(1)
    HeatMap(R)
    %hold on
    %figure(2)
    %plot(1:sideLength, R(midPoint,:))
end

%subplot(ceil((length(matchesTimesToPrint)+2)/4),4,2)
%set(gca,'XLim',[1 tend])
%plot(1:dt:tend, Rstore, 'Color', cmap(find(matchesTimesToPrint),:));
%Xlabel('t')
%ylabel('R(source)')
%hold on

% This is a while loop useful for determining whether the system has
% converged to a steady state yet.
if(0)
    notConverged = 1;
    R = zeros(sideLength,sideLength);
    numIterations = 1;
    while(1)
        Rnext = R + dt.*(DS*([R(:,2:end) R(:,1)] + [R(:,end) R(:,1:end-1)] - ...
            4*R + [R(2:end, :); R(1,:)] + [R(end,:); R(1:end-1, :)])...
            /(dx)^2 + kprod*X);
        if(max(max(abs(R-Rnext))) < .0001*max(max(abs(R))))
            notConverged = 0;
        end
        numIterations = numIterations+1;
        if(mod(numIterations,1000) == 0)
            numIterations
            abs(R(midPoint,midPoint) - R(sideLength,sideLength))
            R(midPoint,midPoint)
            R(sideLength,sideLength)
        end
        R = Rnext;
    end
end