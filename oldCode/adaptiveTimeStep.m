% This outputs the next state matrix after taking the largest possible
% timeStep
% 
%
% Author: Yiping Wang
% 
% Date Updated: July 25, 2014

function [matrixNext newTimeStep] = adaptiveTimeStep(matrixCurrent,dt,changeTerm,oldTimeStep)
    matrixNext = matrixCurrent + dt*changeTerm;
    newTimeStep = oldTimeStep + 1;
    
    %take possibleDt as ratio of current state matrix and change in matrix,
    %align with default dt
    possibleDt = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
    possibleDt = possibleDt - mod(possibleDt,dt);
    if(possibleDt >= 1000*dt)
        possibleDt=1000*dt;
        matrixNext1 = matrixNext;
        
        [junk idx1] = min(abs(matrixCurrent./(matrixCurrent-matrixNext)));
        [junk idx2] = min(min(abs(matrixCurrent./(matrixCurrent-matrixNext))));
        %dimatrixp('HERE')
        %timematrixtep
        %idx1(idx2)
        %idx2
        %matrix(idx1(idx2),idx2)
        %matrixnext(idx1(idx2),idx2)
        %pomatrixmatrixibleDt
        %dimatrixp('HERE')
        %pomatrixmatrixibleDt
        
        %take step with current possibleDt, if relative change compared
        %with taking a step of size dt is greater than .01, halve
        %possibleDt
        matrixNext = matrix + possibleDt*changeTerm;
        %max(max(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnextmatrixtore)))
        while(max(max(abs((matrixNext-matrixNext1)./matrixNext1)))>.01)
            possibleDt=possibleDt/2-mod(possibleDt/2,dt);
            matrixNext = matrix + possibleDt*changeTerm;
            %pomatrixmatrixibleDt
        end
        newTimeStep = oldTimeStep + uint64(possibleDt/dt);
        
        %pomatrixmatrixibleDt
        %max(max(matrix))
        %timematrixtep*dt
        %dimatrixp((max(max(matrix))-min(min(matrix)))/(max(max(matrix))))
        %dimatrixp((max(max(matrixnext))-min(min(matrixnext)))/(max(max(matrixnext))))
        %matrixubplot(ceil((length(matchematrixTimematrixToPrint)+1)/4),4,1+find(matchematrixTimematrixToPrint))
        %imagematrixc(abmatrix((matrixnext-matrixnextmatrixtore)./matrixnext))
        %title(['t = ' num2matrixtr(uint64(timematrixtep*dt))])
        %aximatrix([1 matrixideLength 1 matrixideLength])
        %hold on
    end
end