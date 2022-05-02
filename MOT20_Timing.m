% JLJ
% This script is to time all trackers by running them 10 times on each
% video.

% -------------------------------------------------------------------------
% User settings

clc
clear all
close all

% Get filenames from folder that we will analyze
Directory = 'MOTChallengeEvalKit-master/data/MOT20Labels/train';
Files = dir('MOTChallengeEvalKit-master/data/MOT20Labels/train');
% Files = Files(3:length(Files)) % Get rid of trash entries - Windows version
Files = Files(3:length(Files)) % Get rid of trash entries - Mac version
% -------------------------------------------------------------------------
nTrials=11; % Ignore the first run for each tracker - sometimes it is a lot slower than the rest
% -------------------------------------------------------------------------
%% Tracker Settings
% Flags for which trackers to run
RunOpts.CTv0 = true;
RunOpts.CTv1 = false;
RunOpts.IOU = true;
RunOpts.SORT = true;
RunOpts.OCSORT = true;

DetConf_Thresh = -0.1; % Use all detections for all trackers

% CTv0
CT_Settings.BinWidth = 10; % Could be increased in multiples of 10
CT_Settings.minpts = 2; % DBSCAN parameter, DONT CHANGE
CT_Settings.eps = 16; 
CT_Settings.TrackletMerge_OverlapRatio = 0.4; % "IOU Thresh" above
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


% CTv1 Track Merging Parameters
CT_Settings.Horizon = 60;
CT_Settings.TrackMerge_OverlapRatio = 0;
CT_Settings.VelocityHist = 30;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Read detections from file 
% Det = [f, id, TLx, TLy, W, H, Conf]
for u=1:nTrials
    fprintf('Test: %s\n',num2str(u))
    for i=1:length(Files)
        % -------------------------------------------------------------------------
        name{i} = Files(i).name;
        disp(name{i})
        InputVidName = sprintf('%s.mp4', name{i}(1:8));
        Det = readmatrix(sprintf('%s/%s/det/det.txt', Directory, name{i}));
        Det = sortrows(Det, 1);
        nFrames = max(Det(:,1));

    %     Det = Det(Det(:,7) >  DetConf_Thresh, :);
        Det(:,7) = 1; % MOT20 Provided Detections are only 1 or 0 - 
    %               a confidence of 0 can mess things up so to keep things consisent all detections have conf=1
        % Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz]
        FV = [Det(:,3:4)+Det(:,5:6)/2, ...
                        Det(:,1), Det(:,3:6), ...
                        Det(:,5).*Det(:,6)];
        % ---------------------------------------------------------------------
        % CTv0
        if RunOpts.CTv0 == true
            [~, ~, ~, ~, Time_CTv0(u,i)] = CT_Tracker_v0s(FV, nFrames, 0, CT_Settings);
%             [~, ~, ~, ~, Time_CTv0p(u,i)] = CT_Tracker_v0p(FV, nFrames, 0, CT_Settings);
        end
        % ---------------------------------------------------------------------
        % CTv1
        if RunOpts.CTv1 == true
            [~, ~, ~, ~, Time_CTv1(u,i)] = CT_Tracker_v1s(FV, nFrames, 0, CT_Settings);
        end
        % ---------------------------------------------------------------------
        % IOU Tracker
        if RunOpts.IOU == true
            [~, Time_IOU(u,i)] = run_IOU(Det, IOU_Settings);
        end
        % ---------------------------------------------------------------------
        % SORT Tracker
        if RunOpts.SORT == true
            [~, Time_SORT(u,i)] = run_SORT(Det, SORT_Settings);
        end
        % ---------------------------------------------------------------------
        % OCSORT Tracker
        if RunOpts.OCSORT == true
            [~, Time_OCSORT(u,i)] = run_OCSORT(Det, OCSORT_Settings);
        end
        % -------------------------------------------------------------------------
    end
end

% Get filenames from folder that we will analyze
Directory = 'MOTChallengeEvalKit-master/data/MOT20Labels/test';
Files = dir('MOTChallengeEvalKit-master/data/MOT20Labels/test');
Files = Files(3:length(Files)) % Get rid of trash entries - Windows version
% Files = Files(4:length(Files)) % Get rid of trash entries - Mac version
% -------------------------------------------------------

%% Read detections from file 
% Det = [f, id, TLx, TLy, W, H, Conf]
for u=1:nTrials
    fprintf('Test: %s\n',num2str(u))
    for i=1:length(Files)
        % -------------------------------------------------------------------------
        name{i} = Files(i).name;
        disp(name{i})
        InputVidName = sprintf('%s.mp4', name{i}(1:8));
        Det = readmatrix(sprintf('%s/%s/det/det.txt', Directory, name{i}));
        Det = sortrows(Det, 1);
        nFrames = max(Det(:,1));

    %     Det = Det(Det(:,7) >  DetConf_Thresh, :);
        Det(:,7) = 1; % MOT20 Provided Detections are only 1 or 0 - 
    %               a confidence of 0 can mess things up so to keep things consisent all detections have conf=1
        % Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz]
        FV = [Det(:,3:4)+Det(:,5:6)/2, ...
                        Det(:,1), Det(:,3:6), ...
                        Det(:,5).*Det(:,6)];
        % ---------------------------------------------------------------------
        % ---------------------------------------------------------------------
        % CTv0
        if RunOpts.CTv0 == true
            [~, ~, ~, ~, Time_CTv0(u,i+4)] = CT_Tracker_v0s(FV, nFrames, 0, CT_Settings);
%             [~, ~, ~, ~, Time_CTv0p(u,i)] = CT_Tracker_v0p(FV, nFrames, 0, CT_Settings);
        end
        % ---------------------------------------------------------------------
        % CTv1
        if RunOpts.CTv1 == true
            [~, ~, ~, ~, Time_CTv1(u,i+4)] = CT_Tracker_v1s(FV, nFrames, 0, CT_Settings);
        end
        % ---------------------------------------------------------------------
        % IOU Tracker
        if RunOpts.IOU == true
            [~, Time_IOU(u,i+4)] = run_IOU(Det, IOU_Settings);
        end
        % ---------------------------------------------------------------------
        % SORT Tracker
        if RunOpts.SORT == true
            [~, Time_SORT(u,i+4)] = run_SORT(Det, SORT_Settings);
        end
        % ---------------------------------------------------------------------
        % OCSORT Tracker
        if RunOpts.OCSORT == true
            [~, Time_OCSORT(u,i+4)] = run_OCSORT(Det, OCSORT_Settings);
        end
        % -------------------------------------------------------------------------
    end
end

% Ignore the first run for each tracker - sometimes it is a lot slower than the rest
% writematrix(Time_CTv1,'MOT20_CTv1_Timing')
writematrix(Time_CTv0,'MOT20_CTv0_Timing')
writematrix(Time_IOU,'MOT20_IOU_Timing')
writematrix(Time_SORT,'MOT20_SORT_Timing')
writematrix(Time_OCSORT,'MOT20_OCSORT_Timing')
