% Joachim Lohn-Jaramillo
% Last Update: 4/19/22
% Dartmouth College, PhD Thesis Work
% -------------------------------------------------------------------------
% Our tracker is explained in the following paper, please cite this if you
% use our work.



% This version of CTv0 could be changed to run online with a latency that
% depends on the BinWidth parameter, referred to as tau in the paper. When
% tau=10, we look at 10-frame sequences at a time, which overlap... so
% the 1st sequence contains frames [i through i+10] and the 2nd sequence
% contains frames [i+5 through i+15]).  In this fashion, the tracker would
% run at a 5-10 frame latency depending on how tracklets are connected.

% Note CTv0 expects detections to come in the FV variable with the
% following format: FV = [Cx Cy f TLx TLy Width Height Size]
% Each Track will be output in its own cell array, where each row has the
% following format: Row = [Cx Cy f TLx TLy Width Height Size ID]
function [Tracks, StartPt, EndPt, StationaryFlag, Time] = CT_v0(FV, nFrames, StartF, CT_Settings)
    if ~isempty(FV)
        %%  Split bag (FV) up into bins, with a time with of "CT_Settings.BinWidth"
        %----------------------------------------------------------------------
        % FV = [Cx Cy f TLx TLy Width Height Size]
        tic
        [Bin, Start, End] = BinDetections(FV, nFrames, StartF, CT_Settings.BinWidth);
        time1 = toc;
        %----------------------------------------------------------------------
        %%  Create tracklets for each bin
        %---------------------------------------------------------------------
        % Now we are going to go through each bin, cluster, eliminate outliers,
        % and make sure that there are not multiple objects in each cluster.
        tic
        NumBins = length(Bin);
        MidP = round(CT_Settings.BinWidth/2); % half of CT_Settings.BinWidth, done here to avoid doing it multiple times
        % Pre-allocating
        C         = cell(1,NumBins);
        bBoxStart = cell(1, NumBins);
        bBoxMid   = cell(1,NumBins);
        cStart    = cell(1, NumBins);
        cMid      = cell(1, NumBins);
        ClusterMask=cell(1,NumBins);
        MatchIDs = cell(1,NumBins-1);
        
        ClustersHaveBeenFound = false; % flag for when we find our first detection (sometimes it takes a while)
        MatchDiameter =  3*CT_Settings.eps;
        % Look at each bin of frames independently to cluster and transform
        % clusters into tracklets
        for i=1:NumBins
            FV = Bin{i}(:,1:8);
            if ~isempty(FV)
                idx = dbscan(double(FV(:,1:3)),CT_Settings.eps,CT_Settings.minpts); % run DBSCAN
                labels = unique(idx);
                % Eliminate outliers - noise points are ignored/untracked
                if ismember(-1,labels)
                    FV(idx == labels(1),:) = [];
                    idx(idx == labels(1)) = [];
                    labels(1) = [];
                end
                % -------------------------------------------------------------
                % Go through each proposed cluster and 
                %   - Combine multiple detections from same frame 
                %   - Mark Clusters that do not span more than one frame.
                %   - Fill in any missing detections using a line of best fit
                NumClusters = length(labels);
                ClusterMask{i} = false(NumClusters,1);
                for j=1:NumClusters
                    Temp = FV(idx==labels(j),:); % Temp = points from cluster(j)
                    UniqueF = unique(Temp(:,3));
                    % -------------------------------------------------
                    % We are going to get rid of clusters that 
                    % do not span a large enough segment of time (more than
                    % "CT_Settings.minpts" frames)
                    if length(UniqueF) < CT_Settings.minpts
                        ClusterMask{i}(j) = true;
                        % -------------------------------------------------------------
                        % Make some fake values here, these clusters will be
                        % deleted anyway...
                        % Start point
                        bBoxStart{i}(j,:) = [0 0 0 0];
                        cStart{i}(j,:) = [0 0];

                        % midpoint
                        bBoxMid{i}(j,:) = [0 0 0 0];
                        cMid{i}(j,:) = [0 0];
                        % -------------------------------------------------------------
                        C{i}{j} = Temp;
                    else
                        % -------------------------------------------------------------
                        % If the number of points in the cluster does not match the
                        % number of frames that should be in a bin, go in and
                        % combine detections from same frame.
                        if length(Temp(:,1))~=CT_Settings.BinWidth
                            AlphaCluster = zeros(length(UniqueF),8);
                            % -------------------------------------------------
                            % Combine detections from same frame
                            % Here we are going to check each frame from the
                            % clustered datapoints and see if there are multiple
                            % detections on the same frame, if so -> these
                            % detections will be grouped into a single larger
                            % detection that encapsulates all bboxes. 
                            for p=1:length(UniqueF)
                                Rows = find(Temp(:,3)==UniqueF(p));
                                if length(Rows)>1
                                    bBoxTLCorner = min(Temp(Rows,4:5));
                                    bBoxSize = max((Temp(Rows,4:5) + Temp(Rows,6:7))) - bBoxTLCorner;
                                    Center = bBoxTLCorner + bBoxSize/2;
                                    Size = bBoxSize(1)*bBoxSize(2);
                                    %----------------------
                                    % Cap max size. keep centered around original
                                    % center
                                    if bBoxSize(1) > CT_Settings.MaxTubeDim || bBoxSize(2)>CT_Settings.MaxTubeDim
                                        bBoxSize = [CT_Settings.MaxTubeDim CT_Settings.MaxTubeDim];
                                        bBoxTLCorner = Center - round(CT_Settings.MaxTubeDim/2);
                                        Size = CT_Settings.MaxTubeDim*CT_Settings.MaxTubeDim;
                                    end
                                    %----------------------
                                    AlphaCluster(p,:) = [round(Center) UniqueF(p) bBoxTLCorner bBoxSize Size];
                                else
                                    AlphaCluster(p,:)  = Temp(Rows,:);
                                end
                            end
                            % -----------------------
                            % For tracklets that are too big, I keep the calculated
                            % center of the bounding box, but change the
                            % width/height and top left corner to match the maximum
                            % dimension variable (CT_Settings.MaxTubeDim)
                            SizeMask = find((AlphaCluster(:,6)>CT_Settings.MaxTubeDim | AlphaCluster(:,7)>CT_Settings.MaxTubeDim));
                            if ~isempty(SizeMask)
                                AlphaCluster(SizeMask,4:5) = AlphaCluster(SizeMask,1:2) - round(CT_Settings.MaxTubeDim/2); % Recalc TopLeft bBox location
                                AlphaCluster(SizeMask,6:7) = repmat([CT_Settings.MaxTubeDim, CT_Settings.MaxTubeDim],length(SizeMask), 1); % Recalc TopLeft bBox location
                                AlphaCluster(SizeMask, end) = repmat(CT_Settings.MaxTubeDim^2,length(SizeMask), 1);
                            end
                            % -----------------------
                            % -------------------------------------------------
                            % Fill in missing points....
                            % For each cluster that has missing points,
                            % fit a line of best fit. For missing points, use the
                            % median sized bounding box to fill.
                            MedianBBox = median(AlphaCluster(:,6:7));
                            MedSize = prod(MedianBBox);
                            Range = End(i)+1-Start(i);
                            z = [Start(i):End(i)]; % all frame locations
                            % ------------------------------------------------
                            [XYF] = CTv0_Func_FillMissingPoints(AlphaCluster(:,1:3),z(~ismember(z,AlphaCluster(:,3))));
                            % ------------------------------------------------
                            % Create Cluster: = [Cx Cy f TLx TLy Width Height Size]
                            Temp = [XYF (XYF(:,1:2)-repmat(MedianBBox/2,[Range,1]))...
                                       repmat(MedianBBox,[Range,1]) repmat(MedSize,[Range,1])];
                            % -------------------------------------------------
                        end
                        % -------------------------------------------------------------
                        % Start point
                        bBoxStart{i}(j,:) = Temp(1, 4:7);
                        cStart{i}(j,:) = [Temp(1,1) Temp(1,2)];

                        % midpoint
                        bBoxMid{i}(j,:) = Temp(MidP, 4:7);
                        cMid{i}(j,:) = [Temp(MidP,1) Temp(MidP,2)];
                        % -------------------------------------------------------------
                        C{i}{j} = Temp;
                    end
                    
                end
                if NumClusters>0
                    C{i}(ClusterMask{i})               = [];
                    bBoxStart{i}(ClusterMask{i}, :)    = [];
                    bBoxMid{i}(ClusterMask{i}, :)      = [];
                    cMid{i}(ClusterMask{i}, :)         = [];
                    cStart{i}(ClusterMask{i}, :)       = [];
                end
            else
                C{i} = [];
                cMid{i} = [];
                cStart{i} = [];
                bBoxMid{i} = [];
                bBoxStart{i} = [];
                ClusterMask{i} = 1;
            end
            % -----------------------------------------------------------------
            % -----------------------------------------------------------------
            % -----------------------------------------------------------------
            % -----------------------------------------------------------------
            if i>1
                NumClusters = length(C{i-1}); % !!!!! need to do something when there are no clusters
                if NumClusters==0
                    MatchIDs{i-1} = [];
                else
                    if isempty(cStart{i})
                        % If you get here, that means there are no tracklets to
                        % connect to in the next bin.
                        Match = false(length(cStart{i-1}(:,1)),1);
                        Ind = zeros(size(Match));
                    else
                        [Ind,d] = dsearchn(cStart{i},cMid{i-1});
                        Match1 = d < MatchDiameter; % This is only to limit the number of potential matches to consider
                        % -------------------------------------------------------------
                        % NOTE -----------------------
                        % Here we will verify if a "Match" has a large enough bounding
                        % box overlap to be considered a legit match.  The idea here is
                        % to prevent two bboxes with large size differences from
                        % matching up....
                        % The Matlab function in the below statement unnecessarily calculates bbox
                        % overlap between all possible combinations of both
                        % lists; however, it must be highly optimized
                        % because it runs just as fast as my hand-written
                        % code.
%                         bBoxOverlap = diag(bboxOverlapRatio(bBoxMid{i-1}, bBoxStart{i}(Ind,:)));

                        % Hand-vectorized calculations for bounding box -
                        % (much slower when it is put in a function for
                        % some reason)
                        overlapWidth = min(bBoxMid{i-1}(:,1)+bBoxMid{i-1}(:,3), bBoxStart{i}(Ind,1)+ bBoxStart{i}(Ind,3)) -...
                                          max(bBoxMid{i-1}(:,1), bBoxStart{i}(Ind,1));
                        overlapHeight =  min(bBoxMid{i-1}(:,2)+bBoxMid{i-1}(:,4), bBoxStart{i}(Ind,2)+bBoxStart{i}(Ind,4))-...
                                      max(bBoxMid{i-1}(:,2), bBoxStart{i}(Ind,2));

                        size_intersection = overlapWidth.*overlapHeight;
                        bBoxOverlap =  size_intersection./...
                                        (bBoxMid{i-1}(:,3).* bBoxMid{i-1}(:,4) + ...
                                        bBoxStart{i}(Ind,3).* bBoxStart{i}(Ind,4) - size_intersection);
                        bBoxOverlap(overlapWidth <= 0 | overlapHeight <= 0,1) = 0;
                        % ---------------------------------------------------------
                        Match2 = bBoxOverlap > CT_Settings.TrackletMerge_OverlapRatio;
                        Match = Match1 & Match2;
                        Ind(~Match) = 0;
                        % ---------------------------------------------------------
                        % This section prevents two tracks to connect to the same
                        % cluster.  whichever match has the higher bbox overlap is kept 
                        [Counts, ~] = histcounts(Ind, [0:(max(Ind)+1)]-0.1);
                        RepeatedRows = find(Counts(2:end)>1); % which indices were repeated
                        if ~isempty(RepeatedRows)
                            for j=1:length(RepeatedRows)
                                 Rows = find(Ind == RepeatedRows(j));
                                 tempOverlap = zeros(size(bBoxOverlap));
                                 tempOverlap(Rows, 1) = bBoxOverlap(Rows,1);
                                 [~, tInd] = max(tempOverlap); % tempOverlap was to keep the indices the same between bBoxOverlap
                                 Ind(setdiff(Rows, tInd)) = 0;
                                 Match(setdiff(Rows, tInd)) = false;
                            end
                        end
                    end
                    MatchIDs{i-1} = [[1:NumClusters]' Ind];
                end
            end
            % -----------------------------------------------------------------
            % -----------------------------------------------------------------
            if i>2
                % Build tracks by linking matches
                % -----------------------------------------
                if ~isempty(MatchIDs{i-2}) && ClustersHaveBeenFound == false
                    ClustersHaveBeenFound = true;
                    temp = MatchIDs{i-2}(MatchIDs{i-2}(:,2)~=0,:); % only start the tracks if there was a match (non-matches have a 0 in the second column)
                    for j=1:length(temp(:,1))
                        TrackIDs{j} = [i-2 temp(j,1); i-1 temp(j,2)]; 
                    end
                    ActiveTs = [1:j];
                    nActiveTracks = length(ActiveTs);
                    nTracks = nActiveTracks;

                    % ------------------------------------------
                    % Start New Tracks
                    if ~isempty(MatchIDs{i-1})
                        [Temp,ia] = setdiff(MatchIDs{i-1}(:,1), MatchIDs{i-2}(:,2));
                        for j=1:length(ia)
                            ActiveTs = [ActiveTs nTracks+j];
                            TrackIDs{ActiveTs(end)} = [i-1 Temp(j)];
                        end
                    else
                        if nTracks==0
                            ClustersHaveBeenFound = false;
                        end
                    end
                elseif isempty(MatchIDs{i-2}) && ClustersHaveBeenFound == true
                    % If you are here, then MatchIDs was empty and all current
                    % active tracks are turned off.
                    ActiveTs = [];
                elseif ~isempty(MatchIDs{i-2}) && ClustersHaveBeenFound == true
                    nActiveTracks = length(ActiveTs);
                    % ------------------------------------------
                    % Continue Old Tracks
                    if nActiveTracks>0
                        DeActivateTrack = logical([]);
                        for j=1:nActiveTracks
                            LastCluster = TrackIDs{ActiveTs(j)}(end, 2);
                            ind = find(MatchIDs{i-2}(:,1)==LastCluster);
                            if MatchIDs{i-2}(ind,2) == 0
                                DeActivateTrack(j) = true; 
                            else
                                DeActivateTrack(j) = false; 
                                TrackIDs{ActiveTs(j)} = [TrackIDs{ActiveTs(j)}; [i-1 MatchIDs{i-2}(ind,2)]];
                            end
                        end
                        ActiveTs(DeActivateTrack)=[];
                    end
                    % ------------------------------------------
                    % Start New Tracks
                    if ~isempty(MatchIDs{i-1})
                        nTracks = length(TrackIDs);
                        [Temp,ia] = setdiff(MatchIDs{i-1}(:,1), MatchIDs{i-2}(:,2));
                        for j=1:length(ia)
                            ActiveTs = [ActiveTs nTracks+j];
                            TrackIDs{ActiveTs(end)} = [i-1 Temp(j)];
                        end
                    end
                end
                % -----------------------------------------
            end
        end

        %% ---------------------------------------------------------------------
        % Fill in Tracks using the ID history that was created with TrackID
        % TrackRow = [Cx Cy f TLx TLy Width Height Size]
        NumTracks       = length(TrackIDs);
        Tracks          = cell(1,NumTracks);
        FilteredTrackMask1 = logical([]);
        for T=1:NumTracks
            % TrackIDs{T} has 2 columns, the first is binID, second is clusterID
            if length(TrackIDs{T}(:,1)) >= CT_Settings.MinBins
                FilteredTrackMask1(T) = false;
                for i=1:length(TrackIDs{T}(:,1))
                    Tracks{T} = [Tracks{T}; C{TrackIDs{T}(i,1)}{TrackIDs{T}(i,2)}];
                end
                % -----------------------------------------------------------------
                % Run a short moving average (5-frames) across x,y to make the
                % trajectory smoother.
                Tracks{T}(:,1:2) = movmean(Tracks{T}(:,1:2), 5, 1);
                % -----------------------------------------------------------------
                % I use 'last' here because i want to use the coordinates from the
                % next cluster and not the previous other wise the bounding boxes
                % can lag behind the target... all bc clusters overlap
                [~,UniqueInd,~] = unique(Tracks{T}(:,3),'last');
                Tracks{T} = Tracks{T}(UniqueInd,:);
            else
                FilteredTrackMask1(T) = true;
                % If this condition is met, this track will be deleted and
                % there is no need to actually build it.
            end
        end
        Tracks(FilteredTrackMask1) = [];  

        time2 = toc;
        Time = time1+ time2;
        %----------------------------------------------------------------------
        % Add TrackIDs starting from 1
        % TrackRow = [Cx Cy f TLx TLy Width Height Size]    
        for T=1:length(Tracks)
            Tracks{T} = [Tracks{T} T*ones(length(Tracks{T}(:,1)),1)];
            StartPt(T,:) = Tracks{T}(1, :); % record starting track F.V.
            EndPt(T,:)   = Tracks{T}(end, :); % record last track F.V.
            w_h = median(Tracks{T}(:,6:7));
            StationaryFlag(T) = bboxOverlapRatio([Tracks{T}(1, 1:2)-w_h/2 w_h],[Tracks{T}(end,1:2)-w_h/2 w_h], 'Min') > 0;
        end
        %----------------------------------------------------------------------
    else
        Tracks = {}; 
        StartPt = double.empty(0,8);
        EndPt = double.empty(0,8);
        Time = 0;
    end
    
end
% -------------------------------------------------------------------------

