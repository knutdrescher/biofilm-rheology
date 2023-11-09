function [SeriesID, cellparameters] = TrackInAllFrames(explabel)
tmpdataname = strsplit(explabel, '\');
dataname = tmpdataname{end};

%% Input manual tracking data
TrackingManualProbe_InOneFolder;

%% Track cells between neiboghor two frames
% set diffusive tracking steps. Normally, 100 is an enough number.
maxstep = 100;
%{
TrackOneFrame tracks cells between frame t and t+1.
'setframe' in TrackOneFrame represents frame t+1.
'setframe = 2:4' means cell tracking in frame 1 to 4.
Basically,
frame 1 is an initial state (No flow)
frame 2 is a deformed state (Under shear flow)
frame 3 is a short relaxed state (No flow)
frame 4 is a long relaxed state (No flow)
%}
if ~isempty(cellpos)
    disp(name);
    tic
    %% Track 1 to 2, 2 to 3
    for setframe = 2:3
        [relationcell3, accuracy] = TrackOneFrame(tracklabels, cellparameters, cellImageSize, cellPixelIdxLists, cellpos, setframe, maxstep);
        AutoTrackID{setframe - 1}    = relationcell3;
        AutoTrackIDAcc{setframe - 1} = accuracy;
    end
    toc
    %% Merge Tracking IDs in all frame. Pick up only 1to2, 2to3, 1to3 trajecotries
    SeriesID = AutoTrackID{1};
    for trackframe = 2:length(AutoTrackID)
        for numID = 1: length(SeriesID(:,1))
            targetID = SeriesID(numID, trackframe);
            idx = find(AutoTrackID{trackframe}(:,1) == targetID);
            if ~isempty(idx)
                SeriesID(numID, trackframe+1) = AutoTrackID{trackframe}(idx,2);
            else
                SeriesID(numID, trackframe+1) = 0;
            end
        end
        if ~isempty(SeriesID)
            idx = find(~SeriesID(:, trackframe+1));
            SeriesID(idx,:) = [];
        end
    end
    save([explabel, '\trackedIDs\all_info.mat']);
else
    SeriesID = [];
    cellparameters = [];
end
end
