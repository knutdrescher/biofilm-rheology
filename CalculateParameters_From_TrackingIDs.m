clear;
close all;

topdatafolder = input('Where is the folder including single-cell parameters calculated by BiofilmQ?');
list_exp = dir(topdatafolder);

%% Figure with callback
sum_ratio12 = [];
sum_ratio13 = [];
sum_ratio23 = [];
perfectsum_ratio12 = [];
perfectsum_ratio13 = [];
perfectsum_ratio23 = [];
GoodBad = [];
if exist('AnalyzedData_pca_singlecell_f12f23f13.mat')
    load('AnalyzedData_pca_singlecell_f12f23f13.mat','GoodBad');
end

if isempty(GoodBad)
    %% Want to decide whether each data is good or not
    biofilm_index_13 = 1;
    biofilm_index_23 = 1;
    biofilm_index_12 = 1;
    
    n_goodbad = 1;
    biofilm_index = 1;
    f_n_exp    = 1;
    if exist('tmp_AnalyzedData_pca_singlecell.mat')
        load('tmp_AnalyzedData_pca_singlecell.mat');
        f_n_exp    = n_exp;
        %n_goodbad = n_goodbad - 1;
    end
    
    f = figure('Position', [100 100 2000 1000]);
    set(f, 'WindowButtonUpFcn', @MouseGoodOrNot);
    for n_exp = f_n_exp : length(list_exp)
        path_exp = [list_exp(n_exp).folder,'\',list_exp(n_exp).name];
        list_datamat = dir([path_exp, '\data\*_data.mat']);
        
        pat = digitsPattern;
        extnum = extract(list_exp(f_n_exp).name, pat);
        flowrate = extnum{end-2};
        duration = extnum{end-1};
        
        %% Only data with over 100 trackings
        path_1to3tracks = [path_exp, '\trackedIDs\all_info.mat'];
        load(path_1to3tracks, 'SeriesID');
        if length(SeriesID) > 100
            
            path_datamat1 = [path_exp, '\data\', list_datamat(1).name];
            path_datamat2 = [path_exp, '\data\', list_datamat(2).name];
            path_datamat3 = [path_exp, '\data\', list_datamat(3).name];
            
            load(path_datamat1, 'globalMeasurements');
            if globalMeasurements.Biofilm_Height > 0
                a = struct2table(globalMeasurements);
                load(path_datamat2, 'globalMeasurements');
                b = struct2table(globalMeasurements);
                load(path_datamat3, 'globalMeasurements');
                c = struct2table(globalMeasurements);
                table_globalMeasurements = [a;b;c];
                
                
                %% Calcu
                param_track_f12 = Calculate_Parameters_SingelTrajectories_f1f2(path_exp);
                param_track_f23 = Calculate_Parameters_SingelTrajectories_f2f3(path_exp, param_track_f12);
                param_track_f13 = Calculate_Parameters_SingelTrajectories_f1f3(path_exp, param_track_f12, param_track_f23);
                
                ratio21_f13 = table2array(param_track_f13);
                ratio21_f12 = table2array(param_track_f12);
                ratio21_f23 = table2array(param_track_f23);
                
                Biofilm_z = a.Biofilm_Height;
                Biofilm_xy = sqrt(a.Biofilm_Width * a.Biofilm_Length)/2;
                aspect1 = Biofilm_z/Biofilm_xy;
                
                %% make data line GoodBad
                GoodBad(n_goodbad).exp     = list_exp(n_exp).name;
                %% Check linearity of dz on z13
                subplot(2,3,1);scatter(param_track_f13.z1, param_track_f13.dz);
                title('z vs dz (frame 1 to 3)');
                grid on;
                set(gca, 'FontSize',15);
                subplot(2,3,4);
                plot([param_track_f13.x1, param_track_f13.x1 + param_track_f13.dx]',...
                    [param_track_f13.z1, param_track_f13.z1 + param_track_f13.dz]');
                title('x vs y (frame 1 to 3)');
                axis equal;
                xlim([-30 30]);
                ylim([0 30]);
                set(gca, 'FontSize',15);
                
                %% Check linearity of dz on z f12
                subplot(2,3,2);scatter(param_track_f12.z1, param_track_f12.dz);
                title('z vs dz (frame 1 to 2)');
                grid on;
                set(gca, 'FontSize',15);
                subplot(2,3,5);
                plot([param_track_f12.x1, param_track_f12.x1 + param_track_f12.dx]',...
                    [param_track_f12.z1, param_track_f12.z1 + param_track_f12.dz]');
                title('x vs y (frame 1 to 2)');
                axis equal;
                xlim([-30 30]);
                ylim([0 30]);
                set(gca, 'FontSize',15);
                
                %% Check linearity of dz on z f23
                subplot(2,3,3);scatter(param_track_f23.z1, param_track_f23.dz);
                title('z vs dz (frame 2 to 3)');
                grid on;
                set(gca, 'FontSize',15);
                subplot(2,3,6);
                plot([param_track_f23.x1, param_track_f23.x1 + param_track_f23.dx]',...
                    [param_track_f23.z1, param_track_f23.z1 + param_track_f23.dz]');
                title('x vs y (frame 2 to 3)');
                axis equal;
                xlim([-30 30]);
                ylim([0 30]);
                set(gca, 'FontSize',15);
                
                %% Decision for 13
                f.UserData = 0;
                disp([num2str(n_exp), '/', num2str(length(list_exp))]);
                disp(list_exp(f_n_exp).name);
                disp(['Flow rate: ', flowrate]);
                disp(['Duration: ', duration]);
                disp('Now decision for frame 1 to 3');
                waitfor(f, 'UserData');
                disp('Decision for frame 1 to 3 is done!');
                if f.UserData == 1
                    % good for f13
                    ratio21_f13(:,1) = ratio21_f13(:,1)/Biofilm_xy;
                    ratio21_f13(:,2) = ratio21_f13(:,2)/Biofilm_xy;
                    ratio21_f13(:,3) = ratio21_f13(:,3)/Biofilm_z;
                    ratio21_f13(:,4) = ratio21_f13(:,4)/Biofilm_z;
                    
                    array_aspect1 = ones(size(ratio21_f13,1), 1) * aspect1;
                    array_Biofilm_Volume = ones(size(ratio21_f13,1), 1) * a.Biofilm_Volume;
                    array_flowrate = ones(size(ratio21_f13,1), 1) * str2num(flowrate);
                    array_duration = ones(size(ratio21_f13,1), 1) * str2num(duration);
                    array_Biofilm_xy = ones(size(ratio21_f13,1), 1) * Biofilm_xy;
                    array_Biofilm_z = ones(size(ratio21_f13,1), 1) * Biofilm_z;
                    
                    array_biofilm_index13 = ones(length(ratio21_f13), 1) * biofilm_index_13;
                    
                    ratio21_f13 = [ratio21_f13, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                    sum_ratio13 = [sum_ratio13; [ratio21_f13, array_biofilm_index13]];
                    
                    biofilm_index_13 = biofilm_index_13 + 1;
                    
                end
                GoodBad(n_goodbad).goodbad13 = f.UserData;
                
                %% Decision for 12
                f.UserData = 0;
                disp([num2str(n_exp), '/', num2str(length(list_exp))]);
                disp(list_exp(f_n_exp).name);
                disp(['Flow rate: ', flowrate]);
                disp(['Duration: ', duration]);
                disp('Now decision for frame 1 to 2');
                waitfor(f, 'UserData');
                disp('Decision for frame 1 to 2 is done!');
                if f.UserData == 1
                    ratio21_f12(:,1) = ratio21_f12(:,1)/Biofilm_xy;
                    ratio21_f12(:,2) = ratio21_f12(:,2)/Biofilm_xy;
                    ratio21_f12(:,3) = ratio21_f12(:,3)/Biofilm_z;
                    ratio21_f12(:,4) = ratio21_f12(:,4)/Biofilm_z;
                    
                    array_aspect1 = ones(length(ratio21_f12), 1) * aspect1;
                    array_Biofilm_Volume = ones(length(ratio21_f12), 1) * a.Biofilm_Volume;
                    array_flowrate = ones(length(ratio21_f12), 1) * str2num(flowrate);
                    array_duration = ones(length(ratio21_f12), 1) * str2num(duration);
                    array_Biofilm_xy = ones(size(ratio21_f12,1), 1) * Biofilm_xy;
                    array_Biofilm_z = ones(size(ratio21_f12,1), 1) * Biofilm_z;
                    
                    array_biofilm_index12 = ones(length(ratio21_f12), 1) * biofilm_index_12;
                    
                    ratio21_f12 = [ratio21_f12, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                    sum_ratio12 = [sum_ratio12; [ratio21_f12, array_biofilm_index12]];
                    
                    biofilm_index_12 = biofilm_index_12 + 1;
                end
                GoodBad(n_goodbad).goodbad12 = f.UserData;
                
                %% frame 2 to 3
                f.UserData = 0;
                disp([num2str(n_exp), '/', num2str(length(list_exp))]);
                disp(list_exp(f_n_exp).name);
                disp(['Flow rate: ', flowrate]);
                disp(['Duration: ', duration]);
                disp('Now decision for frame 2 to 3');
                waitfor(f, 'UserData');
                disp('Decision for frame 2 to 3 is done!');
                if f.UserData == 1
                    ratio21_f23(:,1) = ratio21_f23(:,1)/Biofilm_xy;
                    ratio21_f23(:,2) = ratio21_f23(:,2)/Biofilm_xy;
                    ratio21_f23(:,3) = ratio21_f23(:,3)/Biofilm_z;
                    ratio21_f23(:,4) = ratio21_f23(:,4)/Biofilm_z;
                    
                    array_aspect1 =    ones(length(ratio21_f23), 1) * aspect1;
                    array_Biofilm_Volume = ones(length(ratio21_f23), 1) * a.Biofilm_Volume;
                    array_flowrate =   ones(length(ratio21_f23), 1) * str2num(flowrate);
                    array_duration =   ones(length(ratio21_f23), 1) * str2num(duration);
                    array_Biofilm_xy = ones(size(ratio21_f23,1), 1) * Biofilm_xy;
                    array_Biofilm_z =  ones(size(ratio21_f23,1), 1) * Biofilm_z;
                    
                    array_biofilm_index23 = ones(length(ratio21_f23), 1) * biofilm_index_23;
                    
                    ratio21_f23 = [ratio21_f23, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                    sum_ratio23 = [sum_ratio23; [ratio21_f23, array_biofilm_index23]];
                    
                    biofilm_index_23 = biofilm_index_23 + 1;
                end
                GoodBad(n_goodbad).goodbad23 = f.UserData;
                
                
                %% Choose only perfect tracking data
                if GoodBad(n_goodbad).goodbad13 == 1 ...
                        && GoodBad(n_goodbad).goodbad12 == 1 ...
                        && GoodBad(n_goodbad).goodbad23 == 1
                    array_biofilm_index12 = ones(length(ratio21_f12), 1) * biofilm_index;
                    array_biofilm_index13 = ones(length(ratio21_f13), 1) * biofilm_index;
                    array_biofilm_index23 = ones(length(ratio21_f23), 1) * biofilm_index;
                    
                    perfectsum_ratio12  = [perfectsum_ratio12; [ratio21_f12, array_biofilm_index12]];
                    perfectsum_ratio13  = [perfectsum_ratio13; [ratio21_f13, array_biofilm_index13]];
                    perfectsum_ratio23  = [perfectsum_ratio23; [ratio21_f23, array_biofilm_index23]];
                    
                    biofilm_index = biofilm_index + 1;
                    
                end
                n_goodbad = n_goodbad + 1;
            end
        end
        
        sum_ratio12 = unique(sum_ratio12,'rows');
        sum_ratio13 = unique(sum_ratio13,'rows');
        sum_ratio23 = unique(sum_ratio23,'rows');
        
        save('tmp_AnalyzedData_pca_singlecell.mat', 'sum_ratio12','sum_ratio23', 'sum_ratio13', ...
            'perfectsum_ratio12', 'perfectsum_ratio13', 'perfectsum_ratio23', ...
            'GoodBad', 'n_exp', 'n_goodbad', ...
            'biofilm_index', 'biofilm_index_12', 'biofilm_index_23', 'biofilm_index_13');
    end
    close all;
    
    sum_ratio12 = unique(sum_ratio12,'rows');
    sum_ratio13 = unique(sum_ratio13,'rows');
    sum_ratio23 = unique(sum_ratio23,'rows');
    perfectsum_ratio12 = unique(perfectsum_ratio12,'rows');
    perfectsum_ratio13 = unique(perfectsum_ratio13,'rows');
    perfectsum_ratio23 = unique(perfectsum_ratio23,'rows');
    
    save('AnalyzedData_pca_singlecell_f12f23f13.mat');
    
else
    %% Already decided whether each data is good or not. Just want to change included parameters
    biofilm_index = 1;
    biofilm_index_13 = 1;
    biofilm_index_23 = 1;
    biofilm_index_12 = 1;
    
    for n_exp = 1 : length(list_exp)
        pat = digitsPattern;
        extnum = extract(list_exp(n_exp).name, pat);
        flowrate = extnum{end-2};
        duration = extnum{end-1};
        disp(n_exp);        
        
        idx_GoodBad = find(strcmp({GoodBad.exp}, list_exp(n_exp).name));
        
        if ~isempty(idx_GoodBad)
            %% Deleted Multiplied slcies or not
            path_exp = [list_exp(n_exp).folder,'\',list_exp(n_exp).name];
            list_datamat = dir([path_exp, '\data\*_data.mat']);
            %% Only data with over 100 trackings
            path_1to3tracks = [path_exp, '\trackedIDs\all_info.mat'];
            
            if exist(path_1to3tracks)
                load(path_1to3tracks, 'SeriesID');
                if length(SeriesID) > 100
                    
                    path_datamat1 = [path_exp, '\data\', list_datamat(1).name];
                    path_datamat2 = [path_exp, '\data\', list_datamat(2).name];
                    path_datamat3 = [path_exp, '\data\', list_datamat(3).name];
                    
                    load(path_datamat1, 'globalMeasurements','goodObjects');
                    if globalMeasurements.Biofilm_Height > 0
                        a = struct2table(globalMeasurements);
                        load(path_datamat2, 'globalMeasurements');
                        b = struct2table(globalMeasurements);
                        load(path_datamat3, 'globalMeasurements');
                        c = struct2table(globalMeasurements);
                        table_globalMeasurements = [a;b;c];
                        
                        %% Calcu
                        param_track_f12 = Calculate_Parameters_SingelTrajectories_f1f2(path_exp);
                        param_track_f23 = Calculate_Parameters_SingelTrajectories_f2f3(path_exp, param_track_f12);
                        param_track_f13 = Calculate_Parameters_SingelTrajectories_f1f3(path_exp, param_track_f12, param_track_f23);
                        
                        ratio21_f13 = table2array(param_track_f13);
                        ratio21_f12 = table2array(param_track_f12);
                        ratio21_f23 = table2array(param_track_f23);
                        
                        Biofilm_z = a.Biofilm_Height;
                        Biofilm_xy = sqrt(a.Biofilm_Width * a.Biofilm_Length)/2;
                        aspect1 = Biofilm_z/Biofilm_xy;
                        %% Decision for 13
                        if GoodBad(idx_GoodBad).goodbad13 == 1
                            % good for f13
                            ratio21_f13(:,1) = ratio21_f13(:,1)/Biofilm_xy;
                            ratio21_f13(:,2) = ratio21_f13(:,2)/Biofilm_xy;
                            ratio21_f13(:,3) = ratio21_f13(:,3)/Biofilm_z;
                            ratio21_f13(:,4) = ratio21_f13(:,4)/Biofilm_z;
                            
                            array_aspect1 = ones(size(ratio21_f13,1), 1) * aspect1;
                            array_Biofilm_Volume = ones(size(ratio21_f13,1), 1) * a.Biofilm_Volume;
                            array_flowrate = ones(size(ratio21_f13,1), 1) * str2num(flowrate);
                            array_duration = ones(size(ratio21_f13,1), 1) * str2num(duration);
                            array_Biofilm_xy = ones(size(ratio21_f13,1), 1) * Biofilm_xy;
                            array_Biofilm_z = ones(size(ratio21_f13,1), 1) * Biofilm_z;
                            
                            array_biofilm_index13 = ones(length(ratio21_f13), 1) * biofilm_index_13;
                            
                            ratio21_f13 = [ratio21_f13, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                            sum_ratio13 = [sum_ratio13; [ratio21_f13, array_biofilm_index13]];
                            
                            biofilm_index_13 = biofilm_index_13 + 1;
                            
                        end
                        %% Decision for 12
                        if GoodBad(idx_GoodBad).goodbad12 == 1
                            ratio21_f12(:,1) = ratio21_f12(:,1)/Biofilm_xy;
                            ratio21_f12(:,2) = ratio21_f12(:,2)/Biofilm_xy;
                            ratio21_f12(:,3) = ratio21_f12(:,3)/Biofilm_z;
                            ratio21_f12(:,4) = ratio21_f12(:,4)/Biofilm_z;
                            
                            array_aspect1 = ones(length(ratio21_f12), 1) * aspect1;
                            array_Biofilm_Volume = ones(length(ratio21_f12), 1) * a.Biofilm_Volume;
                            array_flowrate = ones(length(ratio21_f12), 1) * str2num(flowrate);
                            array_duration = ones(length(ratio21_f12), 1) * str2num(duration);
                            array_Biofilm_xy = ones(size(ratio21_f12,1), 1) * Biofilm_xy;
                            array_Biofilm_z = ones(size(ratio21_f12,1), 1) * Biofilm_z;
                            
                            array_biofilm_index12 = ones(length(ratio21_f12), 1) * biofilm_index_12;
                            
                            ratio21_f12 = [ratio21_f12, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                            sum_ratio12 = [sum_ratio12; [ratio21_f12, array_biofilm_index12]];
                            
                            biofilm_index_12 = biofilm_index_12 + 1;
                            
                        end
                        
                        %% frame 2 to 3
                        if GoodBad(idx_GoodBad).goodbad23 == 1
                            ratio21_f23(:,1) = ratio21_f23(:,1)/Biofilm_xy;
                            ratio21_f23(:,2) = ratio21_f23(:,2)/Biofilm_xy;
                            ratio21_f23(:,3) = ratio21_f23(:,3)/Biofilm_z;
                            ratio21_f23(:,4) = ratio21_f23(:,4)/Biofilm_z;
                            
                            array_aspect1 =    ones(length(ratio21_f23), 1) * aspect1;
                            array_Biofilm_Volume = ones(length(ratio21_f23), 1) * a.Biofilm_Volume;
                            array_flowrate =   ones(length(ratio21_f23), 1) * str2num(flowrate);
                            array_duration =   ones(length(ratio21_f23), 1) * str2num(duration);
                            array_Biofilm_xy = ones(size(ratio21_f23,1), 1) * Biofilm_xy;
                            array_Biofilm_z =  ones(size(ratio21_f23,1), 1) * Biofilm_z;
                            
                            array_biofilm_index23 = ones(length(ratio21_f23), 1) * biofilm_index_23;
                            
                            ratio21_f23 = [ratio21_f23, array_Biofilm_xy, array_Biofilm_z, array_aspect1, array_Biofilm_Volume, array_flowrate, array_duration];
                            sum_ratio23 = [sum_ratio23; [ratio21_f23, array_biofilm_index23]];
                            
                            biofilm_index_23 = biofilm_index_23 + 1;
                        end
                        
                        
                        %% Choose only perfect tracking data
                        if GoodBad(idx_GoodBad).goodbad13 == 1 ...
                                && GoodBad(idx_GoodBad).goodbad12 == 1 ...
                                && GoodBad(idx_GoodBad).goodbad23 == 1
                            array_biofilm_index12 = ones(length(ratio21_f12), 1) * biofilm_index;
                            array_biofilm_index13 = ones(length(ratio21_f13), 1) * biofilm_index;
                            array_biofilm_index23 = ones(length(ratio21_f23), 1) * biofilm_index;
                            
                            perfectsum_ratio12  = [perfectsum_ratio12; [ratio21_f12, array_biofilm_index12]];
                            perfectsum_ratio13  = [perfectsum_ratio13; [ratio21_f13, array_biofilm_index13]];
                            perfectsum_ratio23  = [perfectsum_ratio23; [ratio21_f23, array_biofilm_index23]];
                            
                            
                            bidx_expname(biofilm_index).path_exp = path_exp;
                            bidx_expname(biofilm_index).name_exp = list_exp(n_exp).name;
                            bidx_expname(biofilm_index).n_label1 = length(goodObjects);
                            bidx_expname(biofilm_index).n_tracked = length(ratio21_f12);
                            bidx_expname(biofilm_index).ratio_tracked_label1 = length(ratio21_f12)/length(goodObjects);
                            
                            
                            biofilm_index = biofilm_index + 1;
                            
                        end
                    end
                end
                
                sum_ratio12 = unique(sum_ratio12,'rows');
                sum_ratio13 = unique(sum_ratio13,'rows');
                sum_ratio23 = unique(sum_ratio23,'rows');
                
            end
        end
    end
    sum_ratio12 = unique(sum_ratio12,'rows');
    sum_ratio13 = unique(sum_ratio13,'rows');
    sum_ratio23 = unique(sum_ratio23,'rows');
    perfectsum_ratio12 = unique(perfectsum_ratio12,'rows');
    perfectsum_ratio13 = unique(perfectsum_ratio13,'rows');
    perfectsum_ratio23 = unique(perfectsum_ratio23,'rows');
    
    %% Parameter names
    coloredparams = param_track_f13.Properties.VariableNames;
    
    coloredparams{24} = 'Biofilm_xy Frame1';
    coloredparams{25} = 'Biofilm_z Frame1';
    coloredparams{26} = 'Biofilm Aspect Ratio Frame1';
    coloredparams{27} = 'Biofilm_Volume Frame1';
    coloredparams{28} = 'Flow rate';
    coloredparams{29} = 'Duration';
    coloredparams{30} = 'Biofilm_Index';
    
    save('AnalyzedData_pca_singlecell_f12f23f13.mat');
end


