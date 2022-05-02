function [] = MakeVideo_ForegroundDetections(OriginalVideoName, OutputVideoName, Bag, ScaleF, GMM, ShowBlobsInGreen, gVideoFont, gLineWidth)
    OldVid = VideoReader(OriginalVideoName); % video handle
    detector = vision.ForegroundDetector('NumGaussians',GMM.nGauss,'MinimumBackgroundRatio', GMM.BackRatio,'NumTrainingFrames',GMM.nTrainingFrames,'LearningRate', GMM.LearningRate);
    NewVid = VideoWriter(OutputVideoName,'MPEG-4');
    NewVid.FrameRate = 30;
    open(NewVid);

    i=1;
    while hasFrame(OldVid)
        frame = readFrame(OldVid);
        gF = rgb2gray(frame);
        if ShowBlobsInGreen == true
            % The reason for the if/else statements is to build the bF_Old mat
            if i==1
                [fHeight, fWidth] = size(gF);
                gF_small = gF(1:ScaleF:fHeight, 1:ScaleF:fWidth);
                bF_Old = false(size(gF_small));
                [bF, bF_Old] = SmallMorph(detector(gF_small), bF_Old);
                bF = imresize(bF,[fHeight, fWidth],'nearest');

            else
                tic
                [bF, bF_Old] = SmallMorph(detector(gF(1:ScaleF:fHeight, 1:ScaleF:fWidth)), bF_Old);
                bF = imresize(bF,[fHeight, fWidth],'nearest');
                bF_Old(:,:,3) = [];
            end
            %------------------------------------------------------
            % bF is the binary frame... this stuff is to make that
            % visualization where I show motion in green
            frame(:,:,1) = gF;
            frame(:,:,2) = uint8(gF + uint8(bF*255));
            frame(:,:,3) = gF;
        else
            % The reason for the if/else statements is to build the bF_Old mat
            if i==1
                [fHeight, fWidth] = size(gF);
                bF_Old = false(round(fHeight/ScaleF), round(fWidth/ScaleF));

                gF_small = gF(1:ScaleF:fHeight, 1:ScaleF:fWidth);
                [bF, bF_Old] = SmallMorph(detector(gF_small), bF_Old);
                bF = imresize(bF,[fHeight, fWidth],'nearest');

            else
                [bF, bF_Old] = SmallMorph(detector(gF(1:ScaleF:fHeight, 1:ScaleF:fWidth)), bF_Old);
                bF = imresize(bF,[fHeight, fWidth],'nearest');
                bF_Old(:,:,3) = [];
            end
            frame = uint8(255*bF);
        end
        %------------------------------------------------------
        Rows = find(Bag(:,3) ==i);
        if ~isempty(Rows)
            frame = insertShape(frame, 'Rectangle', Bag(Rows,4:7),'Color','Red','LineWidth',gLineWidth);
        end
        frame = insertText(frame,[1 1],sprintf('Frame %s',num2str(i)),'FontSize', gVideoFont);
        writeVideo(NewVid,frame);  
        i=i+1;
    end
    close(NewVid)
end

