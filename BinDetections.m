% -------------------------------------------------------------------------
% This function splits up the bag of detections for the entire video into Bins
% of width "BinWidth," which has to date often been 10 frames.
% StartF is an offset here in case the FV list didnt start at frame 1.
function [Bin, Start, End] = BinDetections(FV, nFrames, StartF, BinWidth)
    BinOverlap = 0.5;
    Increment = round(BinWidth*BinOverlap);
    NumBins = floor(nFrames/BinWidth/BinOverlap);
    for i=1:NumBins
        % Start and End are beginning and end frames of the Bin
        if i==NumBins
            Start(i) = (i-1)*Increment+1+StartF;
            End(i) = nFrames+StartF;
        elseif i==1
            Start(i) = (i-1)*BinWidth+1+StartF;
            End(i) = i*BinWidth+StartF;
        else
            Start(i) = (i-1)*Increment+1+StartF;
            End(i) = (i-1)*Increment + BinWidth+StartF;
        end
        BinMask = FV(:,3) >=Start(i) & FV(:,3)<=End(i);
        if isempty(BinMask)
            Bin{i} = [];
        else
            Bin{i} = FV(BinMask,:);
        end
    end    
end

