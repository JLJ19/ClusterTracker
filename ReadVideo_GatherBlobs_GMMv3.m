function [FV, nFrames, Time_GMM, Time_BlobExtraction, fWidth, fHeight] = ReadVideo_GatherBlobs_GMMv3(VidName, GMM, BlobFilt, ScaleF)
% This function reads a video, performs background subtraction, and returns
% a matrix populated with bounding boxes.
% EachRow = [Cx Cy f TLx TLy Width Height Size]
    FV = [];
    F = [];
    Time_GMM = 0;
    Time_BlobExtraction = 0;
    
    tic
    detector = vision.ForegroundDetector('NumGaussians',GMM.nGauss,'MinimumBackgroundRatio', GMM.BackRatio,'NumTrainingFrames',GMM.nTrainingFrames, 'LearningRate', GMM.LearningRate);
    Time_GMM = toc;
        
    Vid = VideoReader(VidName); % video handle
    i=1;
    while hasFrame(Vid)
        frame = readFrame(Vid);
        % The reason for the if/else statements is to build the bF_Old mat
        if i==1
            [fHeight, fWidth,~] = size(frame);
            tic
            gF_small = rgb2gray(frame(1:ScaleF:fHeight, 1:ScaleF:fWidth,:));
            bF_Old = false(size(gF_small));
            [bF, bF_Old] = SmallMorph(detector(gF_small), bF_Old);
        else
            tic
            [bF, bF_Old] = SmallMorph(detector(rgb2gray(frame(1:ScaleF:fHeight, 1:ScaleF:fWidth,:))), bF_Old); % faster to sample frame first before doing grayscale conversion
            bF_Old(:,:,3) = [];
        end
        temp1 = toc;
        %------------------------------------------------------
        % Find all blobs in the video and store in a gigantic bag
        tic
        temp_FV = regionprops(bF,'BoundingBox');
        temp_F = i*ones(length(temp_FV),1);
        
        FV = [FV; temp_FV];
        F = [F; temp_F];
        temp2 = toc;
        %------------------------------------------------------
        i=i+1;
        Time_GMM = Time_GMM+temp1;
        Time_BlobExtraction = Time_BlobExtraction+temp2;
    end
    nFrames = i-1;
    
    tic
    if ~isempty(FV)
        FV = struct2cell(FV)';
        FV = cell2mat(FV)*ScaleF;
        FV = round([FV(:,1:2) + FV(:,3:4)/2 F FV prod(FV(:,3:4),2)]);% FV = [Cx Cy f TLx TLy Width Height Size]

        % Top Down Filter !!! (Only motion within size boundaries)    
        Area = FV(:,6).*FV(:,7);
        Mask = (Area < BlobFilt.MaxObjectSize) & (Area > BlobFilt.MinObjectSize);
        FV = FV(Mask, :);

        % Top Down Filter - No Strange Bounding Box Aspect Ratios
        HeightToWidth = FV(:,7)./FV(:,6);
        Mask = (HeightToWidth > BlobFilt.MinAR) & (HeightToWidth < BlobFilt.MaxAR);
        FV = FV(Mask, :);
    else
        FV = double.empty(0, 8);
    end
    temp2 = toc;
    Time_BlobExtraction = Time_BlobExtraction + temp2;
end

