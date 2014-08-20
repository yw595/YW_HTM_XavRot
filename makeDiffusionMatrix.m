function [diffusionMatrix] = makeDiffusionMatrix(matrix,dx,Dmatrix,timeStep,isCarbon)
    midPoint=50;
    diffusionMatrix = zeros(size(matrix,1),size(matrix,2));
    diffusionMatrix(2:end-1,2:end-1) = Dmatrix*(matrix(2:end-1,1:end-2) + matrix(2:end-1,3:end) - ...
        4*matrix(2:end-1,2:end-1) + matrix(1:end-2,2:end-1) + matrix(3:end,2:end-1))/(dx)^2;
    if(timeStep<10 && isCarbon)
        %diffusionMatrix(midPoint-5:midPoint+5,midPoint-5:midPoint+5)
        if(timeStep==2)
            %matrix
        end
    end
    diffusionMatrix(1,2:end-1) = Dmatrix*(matrix(1,1:end-2) + matrix(1,3:end) + matrix(2,2:end-1) - 3*matrix(1,2:end-1))/(dx)^2;
    diffusionMatrix(end,2:end-1) = Dmatrix*(matrix(end,1:end-2) + matrix(end,3:end) + matrix(end-1,2:end-1) - 3*matrix(end,2:end-1))/(dx)^2;
    diffusionMatrix(2:end-1,1) = Dmatrix*(matrix(1:end-2,1) + matrix(3:end,1) + matrix(2:end-1,2) - 3*matrix(2:end-1,1))/(dx)^2;
    diffusionMatrix(2:end-1,end) = Dmatrix*(matrix(1:end-2,end) + matrix(3:end,end) + matrix(2:end-1,end-1) - 3*matrix(2:end-1,end))/(dx)^2;
    diffusionMatrix(1,1) = Dmatrix*(matrix(2,1) + matrix(1,2) - 2*matrix(1,1))/(dx)^2;
    diffusionMatrix(1,end) = Dmatrix*(matrix(2,end) + matrix(1,end-1) - 2*matrix(1,end))/(dx)^2;
    diffusionMatrix(end,1) = Dmatrix*(matrix(end-1,1) + matrix(end,2) - 2*matrix(end,1))/(dx)^2;
    diffusionMatrix(end,end) = Dmatrix*(matrix(end-1,end) + matrix(end,end-1) - 2*matrix(end,end))/(dx)^2;
end