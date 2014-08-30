function [extendedMatrix1 extendedMatrix2] = makeExtendedMatrix(matrix)
    global onPlate horzPlateEdge vertPlateEdge sideLength midPoint pointsBorderToInside1 pointsBorderToInside2
    matrix1=matrix;
%     size(matrix(onPlate & ([zeros(sideLength,1) horzPlateEdge(:,1:end-1)] | ...
%         [horzPlateEdge(:,2:end) zeros(sideLength,1)])))
%     size(matrix1(horzPlateEdge))
    matrix1(horzPlateEdge) = matrix(onPlate & ([zeros(sideLength,1) horzPlateEdge(:,1:end-1)] | ...
        [horzPlateEdge(:,2:end) zeros(sideLength,1)]));
    %matrix1(horzPlatePoint(1,:)) = matrix1(1,midPoint);
    %matrix1(horzPlatePoint(end,:)) = matrix1(end,midPoint);
    matrix2=matrix;
    matrix2(vertPlateEdge) = matrix(onPlate & ([zeros(1,sideLength); vertPlateEdge(1:end-1,:)] | ...
        [vertPlateEdge(2:end,:); zeros(1,sideLength)]));
    %matrix2(vertPlatePoint(:,1)) = matrix2(midPoint,1);
    %matrix2(vertPlatePoint(:,end)) = matrix2(midPoint,end);
    pointsBorder1 = keys(pointsBorderToInside1);
    for i=1:length(pointsBorder1)
        matrix1(pointsBorder1{i}) = matrix1(pointsBorderToInside1(pointsBorder1{i}));
    end
    pointsBorder2 = keys(pointsBorderToInside2);
    for i=1:length(pointsBorder2)
        matrix2(pointsBorder2{i}) = matrix2(pointsBorderToInside2(pointsBorder2{i}));
    end
    extendedMatrix1 = [matrix1(1,1) matrix1(1,:) matrix1(1,end);matrix1(:,1) matrix1 matrix1(:,end); matrix1(end,1) matrix1(end,:) matrix1(end,end)];
    extendedMatrix2 = [matrix2(1,1) matrix2(1,:) matrix2(1,end);matrix2(:,1) matrix2 matrix2(:,end); matrix2(end,1) matrix2(end,:) matrix2(end,end)];
end