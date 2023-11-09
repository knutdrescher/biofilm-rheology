
clear;
%topdatafolder = 'D:\Mechanics_Project\TrackedData_opp_fixedQthre_20220924';
topdatafolder = input('Where is the folder including single-cell parameters calculated by BiofilmQ?');

list_exp = dir([topdatafolder, '\Nz*']);

for n_exp = 1: length(list_exp)
    disp([list_exp(n_exp).name,' is being analyzed!']);
    explabel = [list_exp(n_exp).folder, '\', list_exp(n_exp).name];
    
    if exist([explabel, '\trackedIDs\trackedID.mat']) ~= 2
        if exist([explabel, '\trackedIDs']) ~= 7
            mkdir([explabel, '\trackedIDs']);
        end
        [SeriesID, cellparameters] = TrackInAllFrames(explabel);
        
        if ~isempty(SeriesID)
            disp([list_exp(n_exp).name,' is done!']);
            save([explabel, '\trackedIDs\trackedID.mat'],'SeriesID','cellparameters');
        end
    else
        disp([list_exp(n_exp).name,' is already tracked!']);
    end
end

