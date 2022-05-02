% -------------------------------------------------------------------------
% This function overlays track bounding boxes on top of a Black & White
% Motion video.
function [] = MakeVideo_MotionAndTracks(OriginalVideoName, OutputVideoName, Bag, ScaleF, GMM, BlackWhite,  FrameRange,gVideoFont,gLineWidth)
    OldVid = VideoReader(OriginalVideoName); % video handle
    detector = vision.ForegroundDetector('NumGaussians',GMM.nGauss,'MinimumBackgroundRatio', GMM.BackRatio,'NumTrainingFrames',GMM.nTrainingFrames, 'LearningRate', GMM.LearningRate);
    
    NewVid = VideoWriter(OutputVideoName,'MPEG-4');
    NewVid.FrameRate = 30;
    open(NewVid);
    i=1;
    while hasFrame(OldVid)
        frame = readFrame(OldVid);
        if i==1
            [fHeight, fWidth,~] = size(frame);
        end
        if i>=FrameRange(1)
            if BlackWhite == true
                % The reason for the if/else statements is to build the bF_Old mat
                if i==1
                    [fHeight, fWidth,~] = size(frame);
                    tic
                    gF_small = rgb2gray(frame(1:ScaleF:fHeight, 1:ScaleF:fWidth, :));
                    bF_Old = false(size(gF_small));
                    [bF, bF_Old] = SmallMorph(detector(gF_small), bF_Old);
                    bF = imresize(bF,[fHeight, fWidth],'nearest');
                else
                    tic
                    [bF, bF_Old] = SmallMorph(detector(rgb2gray(frame(1:ScaleF:fHeight, 1:ScaleF:fWidth,:))), bF_Old);
                    bF = imresize(bF,[fHeight, fWidth],'nearest');
                    bF_Old(:,:,3) = [];
                end
                frame = uint8(bF*255);
            end
            %------------------------------------------------------
            Rows = find(Bag(:,3) ==i);
            if ~isempty(Rows)
                frame = insertObjectAnnotation(frame, 'Rectangle', Bag(Rows,4:7), Bag(Rows,9),'Color','Yellow','FontSize', gVideoFont,'LineWidth',gLineWidth);
            end
            frame = insertText(frame, [1, 50], sprintf('Frame: %s', num2str(i)), 'FontSize',gVideoFont);
            writeVideo(NewVid,frame);  
            if i>FrameRange(2)
                break;
            end
        end
        i=i+1
    end
    close(NewVid);
end

