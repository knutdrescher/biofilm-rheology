function T = Calculate_Parameters_SingelTrajectories_f2f3(path, param_track_f12)

path_parametersmat = [path, '\parameters.mat'];
load(path_parametersmat);
scaling_dxy = params.scaling_dxy/1000;

path_1to3tracks = [path, '\trackedIDs\all_info.mat'];
load(path_1to3tracks);

%% Extract single cell parameters
cellpos_1 = reshape([cellparameters{1}.Centroid],  3, [])'*scaling_dxy;
cellpos_2 = reshape([cellparameters{2}.Centroid],  3, [])'*scaling_dxy;
cellpos_3 = reshape([cellparameters{3}.Centroid],  3, [])'*scaling_dxy;

%% Frame 1 params to compare same information as initial 
dist_surface = [cellparameters{1}.Distance_ToSurface_resolution20]';
dist_surface = dist_surface(SeriesID(:,1));
LocalDensity = [cellparameters{1}.Architecture_LocalDensity_range30]';
LocalDensity = LocalDensity(SeriesID(:,1));
LocalNumDensity = [cellparameters{1}.Architecture_LocalNumberDensity_range30]';
LocalNumDensity = LocalNumDensity(SeriesID(:,1));
NematicOrderParameter = [cellparameters{1}.Architecture_NematicOrderParameter_range30]';
NematicOrderParameter = NematicOrderParameter(SeriesID(:,1));

Alignment_Flow = [cellparameters{1}.Alignment_Flow]';
Alignment_Flow = Alignment_Flow(SeriesID(:,1));
Alignment_Radial = [cellparameters{1}.Alignment_Radial]';
Alignment_Radial = Alignment_Radial(SeriesID(:,1));
Alignment_Zaxis = [cellparameters{1}.Alignment_Zaxis]';
Alignment_Zaxis = Alignment_Zaxis(SeriesID(:,1));

%% Frame2
dist_surface_2 = [cellparameters{2}.Distance_ToSurface_resolution20]';
dist_surface_2 = dist_surface_2(SeriesID(:,2));
LocalDensity_2 = [cellparameters{2}.Architecture_LocalDensity_range30]';
LocalDensity_2 = LocalDensity_2(SeriesID(:,2));
LocalNumDensity_2 = [cellparameters{2}.Architecture_LocalNumberDensity_range30]';
LocalNumDensity_2 = LocalNumDensity_2(SeriesID(:,2));
CellDistMin_2 = [cellparameters{2}.Distance_InterCellSpacing_Min_range20]';
CellDistMin_2 = CellDistMin_2(SeriesID(:,2));
CellDistMean_2 = [cellparameters{2}.Distance_InterCellSpacing_Mean_range20]';
CellDistMean_2 = CellDistMean_2(SeriesID(:,2));
NearestDist_2 = [cellparameters{2}.Distance_ToNearestNeighbor]';
NearestDist_2 = NearestDist_2(SeriesID(:,2));
NematicOrderParameter_2 = [cellparameters{2}.Architecture_NematicOrderParameter_range30]';
NematicOrderParameter_2 = NematicOrderParameter_2(SeriesID(:,2));

Alignment_Flow_2 = [cellparameters{2}.Alignment_Flow]';
Alignment_Flow_2 = Alignment_Flow_2(SeriesID(:,2));
Alignment_Radial_2 = [cellparameters{2}.Alignment_Radial]';
Alignment_Radial_2 = Alignment_Radial_2(SeriesID(:,2));
Alignment_Zaxis_2 = [cellparameters{2}.Alignment_Zaxis]';
Alignment_Zaxis_2 = Alignment_Zaxis_2(SeriesID(:,2));

%% Frame3
dist_surface_3 = [cellparameters{3}.Distance_ToSurface_resolution20]';
dist_surface_3 = dist_surface_3(SeriesID(:,3));
LocalDensity_3 = [cellparameters{3}.Architecture_LocalDensity_range30]';
LocalDensity_3 = LocalDensity_3(SeriesID(:,3));
LocalNumDensity_3 = [cellparameters{3}.Architecture_LocalNumberDensity_range30]';
LocalNumDensity_3 = LocalNumDensity_3(SeriesID(:,3));
CellDistMin_3 = [cellparameters{3}.Distance_InterCellSpacing_Min_range20]';
CellDistMin_3 = CellDistMin_3(SeriesID(:,3));
CellDistMean_3 = [cellparameters{3}.Distance_InterCellSpacing_Mean_range20]';
CellDistMean_3 = CellDistMean_3(SeriesID(:,3));
NearestDist_3 = [cellparameters{3}.Distance_ToNearestNeighbor]';
NearestDist_3 = NearestDist_3(SeriesID(:,3));
NematicOrderParameter_3 = [cellparameters{3}.Architecture_NematicOrderParameter_range30]';
NematicOrderParameter_3 = NematicOrderParameter_3(SeriesID(:,3));

Alignment_Flow_3 = [cellparameters{3}.Alignment_Flow]';
Alignment_Flow_3 = Alignment_Flow_3(SeriesID(:,3));
Alignment_Radial_3 = [cellparameters{3}.Alignment_Radial]';
Alignment_Radial_3 = Alignment_Radial_3(SeriesID(:,3));
Alignment_Zaxis_3 = [cellparameters{3}.Alignment_Zaxis]';
Alignment_Zaxis_3 = Alignment_Zaxis_3(SeriesID(:,3));

%% Difference of parameters
d_dist_surface = dist_surface_3 - dist_surface_2;
d_LocalDensity = LocalDensity_3 - LocalDensity_2;
d_LocalNumDensity = LocalNumDensity_3 - LocalNumDensity_2;
d_CellDistMin = CellDistMin_3 - CellDistMin_2;
d_CellDistMean = CellDistMean_3 - CellDistMean_2;
d_NearestDist = NearestDist_3 - NearestDist_2;
d_NematicOrderParameter = NematicOrderParameter_3 - NematicOrderParameter_2;

d_Alignment_Flow = Alignment_Flow_3 - Alignment_Flow_2;
d_Alignment_Radial = Alignment_Radial_3 - Alignment_Radial_2;
d_Alignment_Zaxis = Alignment_Zaxis_3 - Alignment_Zaxis_2;

%% Cell position
x = []; y = []; z = [];
for m = 1: length(SeriesID)
    x(m) = cellpos_1(SeriesID(m,1), 2);
    y(m) = cellpos_1(SeriesID(m,1), 1);
    z(m) = cellpos_1(SeriesID(m,1), 3);
    x23(m, :) = [cellpos_2(SeriesID(m,2), 2), cellpos_3(SeriesID(m,3), 2)];
    y23(m, :) = [cellpos_2(SeriesID(m,2), 1), cellpos_3(SeriesID(m,3), 1)];
    z23(m, :) = [cellpos_2(SeriesID(m,2), 3), cellpos_3(SeriesID(m,3), 3)];
end
%% Extract lower positions in biofilm
n_bottom10per = round(length(SeriesID)/20);

[~, idx] = mink(z23(:,1), n_bottom10per);
devbottom_x = mean(x23(idx,2) - x23(idx,1));
devbottom_y = mean(y23(idx,2) - y23(idx,1));
devbottom_z = mean(z23(idx,2) - z23(idx,1));

x23(:, 2) = x23(:, 2) - devbottom_x;
y23(:, 2) = y23(:, 2) - devbottom_y;
z23(:, 2) = z23(:, 2) - devbottom_z;


%% Calculate displacements
dx = diff(x23,1,2);
dy = diff(y23,1,2);
dz = diff(z23,1,2);

%% Positions
x1 = param_track_f12.x1;
y1 = param_track_f12.y1;
z1 = param_track_f12.z1;

%% Merge all params in a table
T = table(x1, y1, z1, dist_surface, LocalDensity, LocalNumDensity, NematicOrderParameter...
    , Alignment_Flow, Alignment_Radial, Alignment_Zaxis...
    , dx, dy, dz, d_dist_surface, d_LocalDensity, d_LocalNumDensity, d_CellDistMin, d_CellDistMean, d_NearestDist, d_NematicOrderParameter...
    , d_Alignment_Flow, d_Alignment_Radial, d_Alignment_Zaxis);


end


