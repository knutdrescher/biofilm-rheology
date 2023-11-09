function [relationcell3, accuracy] = TrackOneFrame(tracklabels, cellparameters, cellImageSize, cellPixelIdxLists, cellpos, setframe, maxstep)
% This fucntion relates cells in two frames, frame t and frame t+1
% 'setframe' is defined as t+1
% find nearest cell to probe cell in frame 2
tmpcellpos{1} = cellpos;
cuurentcellparameters = struct2table(cellparameters{setframe});
beforecellparameters = struct2table(cellparameters{setframe-1});
cuurentcellCentroid = cuurentcellparameters.Centroid;
beforecellCentroid = beforecellparameters.Centroid;
beforecellPixelIdxLists = cellPixelIdxLists{setframe-1};
k = 1;

while k < maxstep
    countidx = 1;
    
    %% Predict unknown surrounding cell positions from known trajectories
    for tracklabel = 1:length(cellpos(:,1,1))
        devx = cellpos(tracklabel, setframe-1, 1) - cellpos(tracklabel, setframe, 1);
        devy = cellpos(tracklabel, setframe-1, 2) - cellpos(tracklabel, setframe, 2);
        devz = cellpos(tracklabel, setframe-1, 3) - cellpos(tracklabel, setframe, 3);
        tracklabelpos = squeeze(cellpos(tracklabel, setframe, 1:3))';
        dxyz = cuurentcellCentroid - tracklabelpos;
        dist = sum(dxyz.^2, 2);
        
        [tmp idx] = mink(dist,5);
        
        % identical cell
        Inearest(tracklabel,1:4) = [cuurentcellCentroid(idx(1),:), idx(1)];
        % 1st nearest cell
        Fnearest(tracklabel,1:4) = [cuurentcellCentroid(idx(2),:), idx(2)];
        % 2nd nearest cell
        Snearest(tracklabel,1:4) = [cuurentcellCentroid(idx(3),:), idx(3)];
        % 3rd nearest cell
        Tnearest(tracklabel,1:4) = [cuurentcellCentroid(idx(4),:), idx(4)];
        % 4th nearest cell
        Ynearest(tracklabel,1:4) = [cuurentcellCentroid(idx(5),:), idx(5)];
        
        predictedNC(countidx,1) = round(Inearest(tracklabel,1) + devx);
        predictedNC(countidx,2) = round(Inearest(tracklabel,2) + devy);
        predictedNC(countidx,3) = round(Inearest(tracklabel,3) + devz);
        predictedNC(countidx,4) = Inearest(tracklabel,4);
        countidx = countidx + 1;
        
        predictedNC(countidx,1) = round(Fnearest(tracklabel,1) + devx);
        predictedNC(countidx,2) = round(Fnearest(tracklabel,2) + devy);
        predictedNC(countidx,3) = round(Fnearest(tracklabel,3) + devz);
        predictedNC(countidx,4) = Fnearest(tracklabel,4);
        countidx = countidx + 1;
        
        predictedNC(countidx,1) = round(Snearest(tracklabel,1) + devx);
        predictedNC(countidx,2) = round(Snearest(tracklabel,2) + devy);
        predictedNC(countidx,3) = round(Snearest(tracklabel,3) + devz);
        predictedNC(countidx,4) = Snearest(tracklabel,4);
        countidx = countidx + 1;
        
        predictedNC(countidx,1) = round(Tnearest(tracklabel,1) + devx);
        predictedNC(countidx,2) = round(Tnearest(tracklabel,2) + devy);
        predictedNC(countidx,3) = round(Tnearest(tracklabel,3) + devz);
        predictedNC(countidx,4) = Tnearest(tracklabel,4);
        countidx = countidx + 1;
        
        predictedNC(countidx,1) = round(Ynearest(tracklabel,1) + devx);
        predictedNC(countidx,2) = round(Ynearest(tracklabel,2) + devy);
        predictedNC(countidx,3) = round(Ynearest(tracklabel,3) + devz);
        predictedNC(countidx,4) = Ynearest(tracklabel,4);
        countidx = countidx + 1;
    end
    
    %% Find overlap between predicted position and actual cell position
    % look for overlap between predicted cell position and pixels of obeject in previous frame
    % pixel list to xyz coordinates
    relationcell = [];
    for predictedcell = 1:length(predictedNC)
        centx = predictedNC(predictedcell, 1);
        centy = predictedNC(predictedcell, 2);
        centz = predictedNC(predictedcell, 3);
        % Caution!!! Coordinates in X and Y can be exchaned!!!!
        if centx < 1
            centx = 1;
        end
        if centy < 1
            centy = 1;
        end
        if centz < 1
            centz = 1;
        end
        if centy > cellImageSize{setframe-1}(1)
            centy = cellImageSize{setframe-1}(1);
        end
        if centx > cellImageSize{setframe-1}(2)
            centx = cellImageSize{setframe-1}(2);
        end
        if centz > cellImageSize{setframe-1}(3)
            centz = cellImageSize{setframe-1}(3);
        end
        targetidx = sub2ind(cellImageSize{setframe-1}, centy, centx, centz);
        dxyz = beforecellCentroid - [centx, centy, centz];
        %        tmpdist = sqrt(sum(dxyz.^2, 2))/1000*61;  % um
        tmpdist = sum(dxyz.^2, 2);        % pixel^2
        for precellidx = 1: length(beforecellCentroid)
            if tmpdist(precellidx) < 2418.7 % 3 um converted to pixel^2
                overlap = find(beforecellPixelIdxLists{precellidx} == targetidx, 1);
                if ~isempty(overlap)
                    relationcell(predictedcell,1) = precellidx;
                    relationcell(predictedcell,2) = predictedNC(predictedcell, 4);
                    break
                end
            end
            %        preframecell(predictedcell) = 0;
        end
    end
    relationcell(relationcell(:,1) == 0,:) = [];
    
    %% Delete overtracking trajectries
    [C,ia,ic] = unique(relationcell(:,1));
    overlabels = [];
    for overtrack = 1:length(ia)
        tergetn = relationcell(ia(overtrack),1);
        idx = find(relationcell(:,1) == tergetn);
        if length(idx) > 1
            [C2,ia2,ic2] = unique(relationcell(idx,2));
            if length(ia2) > 1
                % If overlaptrackings are different, delete all
                overlabels = [overlabels; idx];
            else
                % If track labels in overlaptrackings are identical, preserve
                % one
                overlabels = [overlabels; idx(2:end)];
            end
        end
    end
    relationcell2 = relationcell;
    relationcell2(overlabels,:) = [];
    [C,ia,ic] = unique(relationcell2(:,2));
    overlabels = [];
    for overtrack = 1:length(ia)
        tergetn = relationcell2(ia(overtrack),2);
        idx = find(relationcell2(:,2) == tergetn);
        if length(idx) > 1
            overlabels = [overlabels; idx];
        end
    end
    relationcell2(overlabels,:) = [];
    relationcell3 = relationcell2;
    
    %% Check status
    Nonvaltrack(k) = length(relationcell(:,1));
    Valtrack(k)    = length(relationcell3(:,1));
    if mod(k, 5) == 0
        disp(['Step number = ',num2str(k),', Non-validateds = ',num2str(Nonvaltrack(k)),', Validateds = ',num2str(Valtrack(k))]);
    end
    repeatnumberNonv = length(find(Nonvaltrack(:) == Nonvaltrack(k)));
    repeatnumberV = length(find(Valtrack(:) == Valtrack(k)));
    
    %% Put tracked labels into variable
    cellpos = tmpcellpos{1};
    for tracklabel = 1:length(relationcell3(:,1))
        for framenumber = 1:2
            focusframe = setframe + framenumber -2;
            cellpos(tracklabel,focusframe, :) = cellparameters{focusframe}(relationcell3(tracklabel, framenumber)).Centroid;
        end
        %    plot3(cellpos2(tracklabel,:, 1),cellpos2(tracklabel,:, 2),cellpos2(tracklabel,:, 3), 'y', 'LineWidth', 2);
    end
    
    %% Calculate accuracy by comparing with manual tracked labels
    valparam(k) = length(cellpos);
    [C,ia,ib] = intersect(tracklabels(:,setframe-1), relationcell3(:,1));
    trackdiff = [tracklabels(ia,setframe), relationcell3(ib,2)];
    [C2,ia,ib] = intersect(tracklabels(ia,setframe), relationcell3(ib,2));
    accuracy(k) = length(C2)/length(C);
    %    disp(['Accuracy: ',num2str(accuracy(k) * 100)]);
    k = k+1;
    tmpcellpos{k} = cellpos;
    
    %% Tracking is done when the tracked cells are saturated
    if repeatnumberNonv > 2 && repeatnumberV > 2
        break
    end
end
end
