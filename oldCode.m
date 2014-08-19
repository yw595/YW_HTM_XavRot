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