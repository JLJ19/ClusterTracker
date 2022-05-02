% JLJ
% This script compares IOU, SORT, OCSORT, and our tracker using background
% subtraction. Put in a video name into the VideoList variable to evaluate.
% For most of the videos from the VIRAT dataset, the included settings for
% the background model, "GMM_Settings" work reasonably well.
% -------------------------------------------------------------------------
% User settings

clc
clear all
close all
u=1;


VideoList = {'VIRAT_S_000203_06_001206_001266'};
% VideoList = {'VIRAT_S_010111_02_000198_000310'};
%% Tracker Settings
% Flags for which trackers to run
RunOpts.CTv0 = true;
RunOpts.IOU = false;
RunOpts.SORT = false;
RunOpts.OCSORT = true;
RunOpts.ApplyFilter = true;
RunOpts.MakeVideos = true;
MinDuration = 60; 
MinDisplacement = 60; 

%% Tracker Settings
% CTv0
CT_Settings.BinWidth = 30; % Changed in multiples of 10
CT_Settings.minpts = 6; 
CT_Settings.eps = 20; 
CT_Settings.TrackletMerge_OverlapRatio = 0.1; % "IOU Thresh" above
CT_Settings.MinBins = 2; 
CT_Settings.MaxTubeDim = 1000; % Max bbox dimension allowed for objects - if this is a huge value it won't change anything

% IOU
IOU_Settings.sigma_l = 0;
IOU_Settings.sigma_h = 0.1;
IOU_Settings.sigma_iou = 0.3;
IOU_Settings.t_min = 1;

% SORT
SORT_Settings.max_age=100;
SORT_Settings.min_hits=1;
SORT_Settings.iou_threshold = 0.3;% iou_threshold = ?? didnt specify on motchallenge.net ... using default 0.3

% OCSORT
% Official results use adaptive thresholding for confidence - here we DO NOT
% var_iou_threshold=0.3, var_delta_t=3, var_inertia=0.2 are the best
% values reported in the paper for OCSORT
OCSORT_Settings.det_thresh=-0.1; 
OCSORT_Settings.max_age=30; % paper did not comment on this value, but 30=default
OCSORT_Settings.min_hits=3; % paper did not comment on this value, but 3=default
OCSORT_Settings.iou_threshold=0.3;
OCSORT_Settings.delta_t=3;
OCSORT_Settings.inertia=0.2;

% -----------------------
% -------------------------------------------------------------------------
% GMM Parameters
ScaleF = 3;        
BlobFilt_Settings.MaxObjectSize = 490000; %700x700 is largest object allowed
BlobFilt_Settings.MinObjectSize = 25; %5x5 is smallest object allowed (typically)
BlobFilt_Settings.MaxAR = 20;
BlobFilt_Settings.MinAR = 1/BlobFilt_Settings.MaxAR;

GMM_Settings.nGauss = 3;
GMM_Settings.BackRatio = 0.6; % usually at 0.6 for VIRAT
GMM_Settings.nTrainingFrames = 40;
GMM_Settings.LearningRate = 0.005; %0.005 is default
% -------------------------------------------------------------------------
for u=1:length(VideoList)
    InputVidType = '.mp4';
    % Visualization Settings
    MakeVideos = false;
    gVideoFont = 40;
    gLineWidth = 5;
    ShowFigs = false;
    BlackWhiteVideo = false; % Show GMM output as green overlay on grayscale video
    MaxVidLength = inf;

    % -------------------------------------------------------------------------
    %% Extract motion blobs
    VideoRootNames = VideoList{u};
    InputVideoName = sprintf('%s%s',VideoRootNames, InputVidType); %may need to change .mp4 file type
    [FV, nFrames, GMM_Time, BlobExtractionTime, fWidth, fHeight] = ReadVideo_GatherBlobs_GMMv3(InputVideoName, GMM_Settings, BlobFilt_Settings, ScaleF);
    if RunOpts.MakeVideos == true
        MakeVideo_ForegroundDetections(InputVideoName, sprintf('OutputVideos/%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), FV, ScaleF, GMM_Settings, false, gVideoFont, gLineWidth)
    end
    Total_BS_Time = (GMM_Time + BlobExtractionTime)
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % CTv0
    if RunOpts.CTv0 == true
        [Tracks, ~, ~, ~, Time_CTv0(u)] = CT_Tracker_v0s(FV, nFrames, 0, CT_Settings);
        Bag = round(cell2mat(Tracks'));
        [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
        if RunOpts.MakeVideos == true
            MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/CTv0-RawTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)
        end
        if RunOpts.ApplyFilter == true
            [Tracks] = FilterBadTracks(Tracks, MinDuration, MinDisplacement);           
            Bag = round(cell2mat(Tracks'));
            [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
            if RunOpts.MakeVideos == true
                MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/CTv0-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)
                % ---------------------------------------------------------------------
            end
        end
    end
    % ---------------------------------------------------------------------
    % IOU Tracker
    if RunOpts.IOU == true
        [Bag, Tracks, Time_IOU(u)] = run_IOU_UsingMyBag(FV, 0, IOU_Settings);
        [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
        if RunOpts.MakeVideos == true
            MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/IOU-RawTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
        end
        if RunOpts.ApplyFilter == true
            [Tracks] = FilterBadTracks(Tracks, MinDuration, MinDisplacement);           
            Bag = round(cell2mat(Tracks'));
            [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
            if RunOpts.MakeVideos == true
                MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/IOU-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
            end
        end
    end
    % ---------------------------------------------------------------------
    % SORT Tracker
    if RunOpts.SORT == true
        [Bag, Tracks, Time_SORT(u)] = run_SORT_UsingMyBag(FV, 0, SORT_Settings);
        [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
        if RunOpts.MakeVideos == true
            MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/SORT-RawTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
        end
        if RunOpts.ApplyFilter == true
            [Tracks] = FilterBadTracks(Tracks, MinDuration, MinDisplacement);    
            Bag = cell2mat(Tracks');
            [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
            if RunOpts.MakeVideos == true
                MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/SORT-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
            end
        end
    end
    % ---------------------------------------------------------------------
    % OCSORT Tracker
    if RunOpts.OCSORT == true
        [Bag, Tracks, Time_OCSORT(u)] = run_OCSORT_UsingMyBag(FV, 0 ,  OCSORT_Settings);
        [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
        if RunOpts.MakeVideos == true
            MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/OCSORT-RawTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
        end
        if RunOpts.ApplyFilter == true
            [Tracks] = FilterBadTracks(Tracks, MinDuration, MinDisplacement);           
            Bag = cell2mat(Tracks');
            [Bag] = FixOutOfBoundCoords_KeepSizeSame(Bag, fHeight, fWidth); 
            if RunOpts.MakeVideos == true
                MakeVideo_MotionAndTracks(InputVideoName, sprintf('OutputVideos/OCSORT-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate), Bag, ScaleF, GMM_Settings, BlackWhiteVideo,  [1 inf] ,gVideoFont,gLineWidth)% ---------------------------------------------------------------------
            end
        end
    end
    % -------------------------------------------------------------------------
    % Make a Combined Video: 3 videos stacked vertically
    VideoCombiner_x3(sprintf('OutputVideos/%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f.mp4',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate),...
        sprintf('OutputVideos/CTv0-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f.mp4',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate),...
        sprintf('OutputVideos/OCSORT-FilteredTracks_%s_GMM%0.0f_Ratio%.2f_ScaleF%0.0f_LR%0.3f.mp4',VideoRootNames, GMM_Settings.nGauss, GMM_Settings.BackRatio, ScaleF, GMM_Settings.LearningRate),...
        sprintf('%s_TrackerComparison.mp4',VideoRootNames),30);
    
end
%%
% -------------------------------------------------------------------------
function [FV] = ConvertMOTCMatrix2MyBag(Det)
    % Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz]
    FV = [Det(:,3:4)+Det(:,5:6)/2, ...
                    Det(:,1), Det(:,3:6), ...
                    Det(:,5).*Det(:,6)];
end

function [Tracks] = FilterBadTracks(Tracks, MinDuration, MinDisplacement)
    [StartPt, EndPt, Stationary] = ExtraTrackInfo(Tracks);
    Duration = EndPt(:,3) - StartPt(:,3);
    TotalDisplacement = sqrt(sum(abs((EndPt(:,1:2) - StartPt(:,1:2))).^2, 2));
    % ---------------------------------------------------------------------
    % TrackRows = [Cx Cy f TLx TLy Width Height Size TubeID TrustedFlag]
    TrustedTracksMask = (Stationary == false) & Duration' > MinDuration &  TotalDisplacement' > MinDisplacement;    
    Tracks = Tracks(TrustedTracksMask);
end

function [StartPt, EndPt, StationaryFlag] = ExtraTrackInfo(Tracks)
    for T=1:length(Tracks)
        StartPt(T,:) = Tracks{T}(1, :); % record starting track F.V.
        EndPt(T,:)   = Tracks{T}(end, :); % record last track F.V.
        w_h = [median(Tracks{T}(:,6)), median(Tracks{T}(:,7))];
        StationaryFlag(T) = bboxOverlapRatio([Tracks{T}(1, 1:2)-w_h/2 w_h],[Tracks{T}(end,1:2)-w_h/2 w_h], 'Min') > 0;
    end
end
