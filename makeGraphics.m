function makeGraphics(arrayAll,timesToPrint,figNum,title,cmap,makeMovie)
    matrixToPlot = zeros(sideLength,sideLength,3);
    set(0, 'CurrentFigure', figure(figNum))
    hold on
    if strcmp(title,'L and R')
        MINVAL=min(min(min(min(arrayAll))));
        MAXVAL=max(max(max(max(arrayAll))));
    elseif strcmp(title,'GFP')
        MINVAL=min(min(min(arrayAll(1,:,:,:))));
        MAXVAL=max(max(max(arrayAll(1,:,:,:))));
    elseif strcmp(title,'X')
        MINVAL=min(min(min(arrayAll(2,:,:,:))));
        MAXVAL=max(max(max(arrayAll(2,:,:,:))));
    else
        MINVAL=min(min(min(arrayAll)));
        MAXVAL=max(max(max(arrayAll)));
    end
    
    for j=1:length(timesToPrint)
        if strcmp(title,'L') || strcmp(title,'R')
            array = arrayAll(:,:,j);
        elseif strcmp(title,'X')
            array = arrayAll(2,:,:,j);
        elseif strcmp(title,'GFP')
            array = arrayAll(1,:,:,j);
            array = array + arrayAll(3,:,:,j)*.1*MAXVAL;
        end
        
        if strcmp(title,'L')
            matrixToPlot(:, :, 1) = uint16(0);
            matrixToPlot(:, :, 2) = uint16(0);
            matrixToPlot(:, :, 3) = array;
        elseif strcmp(title,'R')
            matrixToPlot(:, :, 1) = uint16(0);
            matrixToPlot(:, :, 2) = array;
            matrixToPlot(:, :, 3) = uint16(0);
        elseif strcmp(title,'L and R')
            matrixToPlot(:, :, 1) = uint16(0);
            matrixToPlot(:, :, 2) = arrayAll(1,:,:,j);
            matrixToPlot(:, :, 3) = arrayAll(2,:,:,j);
        elseif strcmp(title,'X')
            matrixToPlot(:, :, 1) = uint16(0);
            matrixToPlot(:, :, 2) = uint16(0);
            matrixToPlot(:, :, 3) = array;
        elseif strcmp(title,'GFP')
            matrixToPlot(:, :, 1) = uint16(0);
            matrixToPlot(:, :, 2) = array;
            matrixToPlot(:, :, 3) = uint16(0);
        end
        
        matrixToPlotScale = uint16( (matrixToPlot - MINVAL) / (MAXVAL - MINVAL) * (2^16-1) );
        if strcmp(title,'X')
            matrixToPlotScale(:,:,1) = arrayAll(1,:,:,j)*(2^16-1);
        end
        if strcmp(title,'GFP')
            matrixToPlotScale(:,:,1) = arrayAll(2,:,:,j).*(~arrayAll(3,:,:,j))*(2^16-1);
        end
        
        subplot(ceil((length(timesToPrint)+1)/4),4,1)
        set(gca,'XLim',[1 sideLength])
        if strcmp(title,'L')
            plot(1:sideLength, matrixToPlot(midPoint,:,3), 'Color', ...
                cmap(j,:));
            xlabel('x')
            ylabel('L')
        elseif strcmp(title,'R')
            plot(1:sideLength, matrixToPlot(midPoint,:,2), 'Color', ...
                cmap(j,:));
            xlabel('x')
            ylabel('R')
        end
        hold on
        
        %plot 2d images of L, R, or both
        if strcmp(title,'L') || strcmp(title,'R')
            subplot(ceil((length(timesToPrint)+1)/4), ...
                4, 1 + j)
        else
            subplot(ceil((length(timesToPrint))/4), ...
                4, j)
        end
        imshow(matrixToPlotScale)
        title(['t = ' num2str(timesToPrint(j))])
        axis([1 sideLength 1 sideLength])
        hold on
        
        if(makeMovie)
            figure(15)
            imshow(matrixToPlotScale)
            drawnow
            if strcmp(title,'GFP')
                print( 15, '-djpeg', sprintf('frames/frmGFP2%d.jpg', j))
            end
            if strcmp(title,'X')
                print( 15, '-djpeg', sprintf('frames/frmX2%d.jpg', j))
            end
        end
    end
    
    if(makeMovie && (strcmp(title,'X') || strcmp(title,'GFP')))
        disp('making ffmpeg movie')
        if strcmp(title,'GFP')
            system('"C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmGFP2%d.jpg frames/GFP2.mov')
        end
        if strcmp(title,'X')
            system('"C:\Users\Yiping Wang\Documents\ffmpeg\bin\ffmpeg" -y -r 10 -i frames/frmX2%d.jpg frames/X2.mov')
        end
        disp('done with ffmpeg movie')
    end
end