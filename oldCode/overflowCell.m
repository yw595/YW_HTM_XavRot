function [LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,colonies,newOverflowCoord1 newOverflowCoord2 newOverflowXVal newOverflowLIKOVal newOverflowRIKOVal ...
    newOverflowLRKOVal newOverflowRRKOVal newOverflowIiVal newOverflowNiVal newOverflowColoniesVal] ... 
    = overflowCell(varargin)

    global onPlate

    %randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
    LIKO = varargin{1};
    RIKO = varargin{2};
    LRKO = varargin{3};
    RRKO = varargin{4};
    Ii = varargin{5};
    Ni = varargin{6};
    X = varargin{7};
    colonies = varargin{8};
    overflowCoord1 = varargin{9};
    overflowCoord2 = varargin{10};
    sideLength = varargin{11};
    timeStep = varargin{12};
    colonyCenterCoords = varargin{13};
    colony = varargin{14};
    first=1;
    if(length(varargin)>14)
        overflowXVal = varargin{15};
        overflowLIKOVal = varargin{16};
        overflowRIKOVal = varargin{17};
        overflowLRKOVal = varargin{18};
        overflowRRKOVal = varargin{19};
        overflowIiVal = varargin{20};
        overflowNiVal = varargin{21};
        overflowColoniesVal = varargin{22};
        first=0;
    end
    if(colony~=0)
        centerCoord1 = colonyCenterCoords(colony,1);
        centerCoord2 = colonyCenterCoords(colony,2);
    end
    
    newOverflowLIKOVal = LIKO(overflowCoord1,overflowCoord2);
    newOverflowRIKOVal = RIKO(overflowCoord1,overflowCoord2);
    newOverflowLRKOVal = LRKO(overflowCoord1,overflowCoord2);
    newOverflowRRKOVal = RRKO(overflowCoord1,overflowCoord2);
    newOverflowColoniesVal = colonies(overflowCoord1,overflowCoord2);
    if(first)
        newOverflowIiVal = Ii(overflowCoord1,overflowCoord2);
        newOverflowNiVal = Ni(overflowCoord1,overflowCoord2);
        newOverflowXVal = X(overflowCoord1,overflowCoord2)/2;
    else
        newOverflowIiVal = Ii(overflowCoord1,overflowCoord2);
        newOverflowNiVal = Ni(overflowCoord1,overflowCoord2);
        newOverflowXVal = X(overflowCoord1,overflowCoord2);
    end
    
    if(1)
    [XDist XNearestBorder] = bwdist(X==0 & onPlate,'chessboard');
    nearestBorderCoord1 = mod(XNearestBorder(overflowCoord1,overflowCoord2)-1,sideLength)+1;
    nearestBorderCoord2 = floor((double(XNearestBorder(overflowCoord1,overflowCoord2)) -1 )/ (sideLength*1.0))+1;
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
        if(newOverflowCoord1>sideLength)
            newOverflowCoord1=1;
        end
        if(newOverflowCoord1<1)
            newOverflowCoord1=sideLength;
        end
        newOverflowCoord2=overflowCoord2+randCoords(numRand,2);
        if(newOverflowCoord2>sideLength)
            newOverflowCoord2=1;
        end
        if(newOverflowCoord2<1)
            newOverflowCoord2=sideLength;
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
        if(newOverflowCoord1>sideLength)
            newOverflowCoord1=1;
        end
        if(newOverflowCoord1<1)
            newOverflowCoord1=sideLength;
        end
        newOverflowCoord2=overflowCoord2+randCoords(numRand,2);
        if(newOverflowCoord2>sideLength)
            newOverflowCoord2=1;
        end
        if(newOverflowCoord2<1)
            newOverflowCoord2=sideLength;
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
        if(newOverflowCoord1>sideLength)
            newOverflowCoord1=1;
        end
        if(newOverflowCoord1<1)
            newOverflowCoord1=sideLength;
        end
        if(newOverflowCoord2>sideLength)
            newOverflowCoord2=1;
        end
        if(newOverflowCoord2<1)
            newOverflowCoord2=sideLength;
        end
    end
    end
    
    if(first)
        Ii(overflowCoord1,overflowCoord2) = Ii(overflowCoord1,overflowCoord2);
        Ni(overflowCoord1,overflowCoord2) = Ni(overflowCoord1,overflowCoord2);
        X(overflowCoord1,overflowCoord2) = X(overflowCoord1,overflowCoord2)/2;
    else
        LIKO(overflowCoord1,overflowCoord2) = overflowLIKOVal;
        RIKO(overflowCoord1,overflowCoord2) = overflowRIKOVal;
        LRKO(overflowCoord1,overflowCoord2) = overflowLRKOVal;
        RRKO(overflowCoord1,overflowCoord2) = overflowRRKOVal;
        colonies(overflowCoord1,overflowCoord2) = overflowColoniesVal;
        Ii(overflowCoord1,overflowCoord2) = overflowIiVal;
        Ni(overflowCoord1,overflowCoord2) = overflowNiVal;
        X(overflowCoord1,overflowCoord2) = overflowXVal;
    end
end