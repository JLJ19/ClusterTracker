function [ ] = VideoCombiner_x3(Name1, Name2, Name3, OutputName,FrameRate)
% JLJ
% This function is meant to open 3 videos and combine them for comparison
% purposes. (stacked vertically)
% Top Left: Vid1
% TR: Vid2
% BL: Vid3
% BR: Vid4
% NOTE: it is ASSUMED all videos have same dimensions
Vid1 = VideoReader(Name1);
Vid2 = VideoReader(Name2);
Vid3 = VideoReader(Name3);
fHeight = Vid1.Height;
fWidth  = Vid1.Width;


NewVid = VideoWriter(OutputName,'MPEG-4');
NewVid.FrameRate = FrameRate;

open(NewVid);
NewFrame(fHeight*2,fWidth*2,3)=uint8(0);
CombinedFrames = 1;
while hasFrame(Vid1)
    f1 = readFrame(Vid1); % open frame from first video
    if hasFrame(Vid2)==1
        f2 = readFrame(Vid2); % open frame from second video
    else
        f2 = zeros(fHeight,fWidth);
    end
    if hasFrame(Vid3)==1
        f3 = readFrame(Vid3); % open frame from third video
    else
        f3 = zeros(fHeight,fWidth);
    end
    NewFrame(1:fHeight, 1:fWidth,:) = f1;
    NewFrame(fHeight+1:fHeight*2, 1:fWidth,:) = f2;
    NewFrame(fHeight*2+1:fHeight*3, 1:fWidth,:) = f3;
    writeVideo(NewVid,NewFrame);
    CombinedFrames = CombinedFrames +1
end
close(NewVid);
end