function [LIKO,RIKO,LRKO,RRKO,Ii,Ni,X,newOverflowCoord1 newOverflowCoord2 newOverflowXVal newOverflowLIKOVal newOverflowRIKOVal ...
    newOverflowLRKOVal newOverflowRRKOVal newOverflowIiVal newOverflowNiVal] ... 
    = overflowCell(varargin)%(overflowCoord1,overflowCoord2,overflowXVal,overflowLIKOVal,overflowRIKOVal,...
    %overflowLRKOVal,overflowRRKOVal,overflowIiVal,overflowNiVal)

    randCoords = [-1 1; -1 0; -1 -1; 0 1; 0 -1; 1 1; 1 0; 1 -1];
    sideLength=100;
    LIKO = varargin{1};
    RIKO = varargin{2};
    LRKO = varargin{3};
    RRKO = varargin{4};
    Ii = varargin{5};
    Ni = varargin{6};
    X = varargin{7};
    overflowCoord1 = varargin{8};
    overflowCoord2 = varargin{9};
    first=1;
    if(length(varargin)>9)
        overflowXVal = varargin{10};
        overflowLIKOVal = varargin{11};
        overflowRIKOVal = varargin{12};
        overflowLRKOVal = varargin{13};
        overflowRRKOVal = varargin{14};
        overflowIiVal = varargin{15};
        overflowNiVal = varargin{16};
        first=0;
    end
    
    newOverflowLIKOVal = LIKO(overflowCoord1,overflowCoord2);
    newOverflowRIKOVal = RIKO(overflowCoord1,overflowCoord2);
    newOverflowLRKOVal = LRKO(overflowCoord1,overflowCoord2);
    newOverflowRRKOVal = RRKO(overflowCoord1,overflowCoord2);
    if(first)
        newOverflowIiVal = Ii(overflowCoord1,overflowCoord2)/2;
        newOverflowNiVal = Ni(overflowCoord1,overflowCoord2)/2;
        newOverflowXVal = X(overflowCoord1,overflowCoord2)/2;
    else
        newOverflowIiVal = Ii(overflowCoord1,overflowCoord2);
        newOverflowNiVal = Ni(overflowCoord1,overflowCoord2);
        newOverflowXVal = X(overflowCoord1,overflowCoord2);
    end
    numRand = randi(8,1);
    newOverflowCoord1=overflowCoord1+randCoords(numRand,1);
    if(newOverflowCoord1>sideLength)
        newOverflowCoord1=1;
    end
    newOverflowCoord2=overflowCoord2+randCoords(numRand,2);
    if(newOverflowCoord2>sideLength)
        newOverflowCoord2=1;
    end
    
    if(first)
        Ii(overflowCoord1,overflowCoord2) = Ii(overflowCoord1,overflowCoord2)/2;
        Ni(overflowCoord1,overflowCoord2) = Ni(overflowCoord1,overflowCoord2)/2;
        X(overflowCoord1,overflowCoord2) = X(overflowCoord1,overflowCoord2)/2;
    else
        LIKO(overflowCoord1,overflowCoord2) = overflowLIKOVal;
        RIKO(overflowCoord1,overflowCoord2) = overflowRIKOVal;
        LRKO(overflowCoord1,overflowCoord2) = overflowLRKOVal;
        RRKO(overflowCoord1,overflowCoord2) = overflowRRKOVal;
        Ii(overflowCoord1,overflowCoord2) = overflowIiVal;
        Ni(overflowCoord1,overflowCoord2) = overflowNiVal;
        X(overflowCoord1,overflowCoord2) = overflowXVal;
    end
end