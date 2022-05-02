% JLJ
% This script evaluates a tracker with a range of parameters.  This is used
% to find the the best set of results

% Assumed Dataset Structure
% MOT20
%       -- train
%               -- MOT20-01
%                       -- det
%                               -- det.txt
%               -- MOT20-02
clc
clear all
close all

% -------------------------------------------------------------------------
% Get filenames from folder that we will analyze
DatasetName = 'MOT20';
Directory = sprintf('Datasets/%s/train', DatasetName);
Files = dir(Directory);
Files = Files(3:length(Files)) % Get rid of non-files, there should be 4 for MOT20 train
% -------------------------------------------------------------------------


%% Tracker Settings
% Flags for which trackers to run
RunOpts.CTv0 = true;
RunOpts.IOU = true;
RunOpts.SORT = true;
RunOpts.OCSORT = true;

DetConf_Thresh = -0.1; % Use all detections for all trackers

% CTv0

CT_Settings.eps = 16; 
CT_Settings.TrackletMerge_OverlapRatio = 0.4; % "IOU Thresh" above
CT_Settings.MinBins = 2;
CT_Settings.BinWidth = 10; % Could be increased in multiples of 10
CT_Settings.MaxTubeDim = 1000; % Max bbox dimension allowed for an object
CT_Settings.minpts = 2; % DBSCAN parameter, DONT CHANGE

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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Read detections from file 
% Det = [f, id, TLx, TLy, W, H, Conf]
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
        [CT_Tracks, ~, ~, ~, Time_CTv0(i)] = CT_Tracker_v0s(FV, nFrames, 0, CT_Settings);
        Bag = round(cell2mat(CT_Tracks'));
        % MyBAG =  [Cx Cy f TLx TLy Width Height TubeID]
        % MOT Data Style = <frame>, <id>, <bb_left>, <bb_top>, <bb_width>, <bb_height>, <conf>, <x>, <y>, <z>
        MotBag_Det = [Bag(:,3), Bag(:,end), Bag(:,4:7), ones(size(Bag(:,1))), -1*ones([length(Bag(:,1)) 3])];
        writematrix(MotBag_Det, sprintf('Output/MOT20-train/CTv0/data/%s.txt', name{i}));
    end
    % ---------------------------------------------------------------------
    % IOU Tracker
    if RunOpts.IOU == true
        [track_result, Time_IOU(i)] = run_IOU(Det, IOU_Settings);
        writematrix(track_result, sprintf('Output/MOT20-train/IOU/data/%s.txt', name{i}));
    end
    % ---------------------------------------------------------------------
    % SORT Tracker
    if RunOpts.SORT == true
        [track_result, Time_SORT(i)] = run_SORT(Det, SORT_Settings);
        writematrix(track_result, sprintf('Output/MOT20-train/SORT/data/%s.txt', name{i}));
    end
    % ---------------------------------------------------------------------
    % OCSORT Tracker
    if RunOpts.OCSORT == true
        [track_result, Time_OCSORT(i)] = run_OCSORT(Det, OCSORT_Settings);
        writematrix(track_result, sprintf('Output/MOT20-train/OCSORT/data/%s.txt', name{i}));
    end
    % -------------------------------------------------------------------------
end
