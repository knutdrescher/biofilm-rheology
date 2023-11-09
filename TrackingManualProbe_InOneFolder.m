%% Load manual tracking data
path_ManualTrackedCellLabelIDs = input('Where is a manually tracked-cell ID?');
load(path_ManualTrackedCellLabelIDs);
%load('ManualTrackedLabel_newOPP.mat');
[~,name,~] = fileparts(explabel);

cellpos = [];
tracklabels = [];
idx_exp = find(strcmp(T.DataName, name));
if ~isempty(idx_exp)
    tracklabels(1) = T.Frame1(idx_exp);
    tracklabels(2) = T.Frame2(idx_exp);
    tracklabels(3) = T.Frame3(idx_exp);
    
    %% load cell parameters
    datapath = [explabel,'\data'];
    datalist = dir([datapath,'\*ch1*.mat']);
    for framenumber = 1: length(datalist)
        framedata = [datalist(framenumber).folder,'\',datalist(framenumber).name];
        load(framedata, 'stats', 'PixelIdxList', 'ImageSize');
        cellparameters{framenumber} = stats;
        cellPixelIdxLists{framenumber} = PixelIdxList;
        cellImageSize{framenumber} = ImageSize;
        clearvars stats PixelIdxList ImageSize
    end
    
    %% make variable which contains cell positions tracked manually
    for tracklabel = 1:length(tracklabels(:,1))
        for framenumber = 1:length(tracklabels(1,:))
            cellpos(tracklabel,framenumber, :) = cellparameters{framenumber}(tracklabels(tracklabel, framenumber)).Centroid;
        end
    end
end

