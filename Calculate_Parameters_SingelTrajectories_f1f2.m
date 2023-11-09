function T = Calculate_Parameters_SingelTrajectories_f1f2(path)

path_parametersmat = [path, '\parameters.mat'];
load(path_parametersmat);
scaling_dxy = params.scaling_dxy/1000;

path_1to3tracks = [path, '\trackedIDs\all_info.mat'];
load(path_1to3tracks);

%% Extract single cell parameters

cellpos_1 = reshape([cellparameters{1}.Centroid],  3, [])'*scaling_dxy;
cellpos_2 = reshape([cellparameters{2}.Centroid],  3, [])'*scaling_dxy;

dist_surface = [cellparameters{1}.Distance_ToSurface_resolution20]';
dist_surface = dist_surface(SeriesID(:,1));
LocalDensity = [cellparameters{1}.Architecture_LocalDensity_range30]';
LocalDensity = LocalDensity(SeriesID(:,1));
LocalNumDensity = [cellparameters{1}.Architecture_LocalNumberDensity_range30]';
LocalNumDensity = LocalNumDensity(SeriesID(:,1));
CellDistMin = [cellparameters{1}.Distance_InterCellSpacing_Min_range20]';
CellDistMin = CellDistMin(SeriesID(:,1));
CellDistMean = [cellparameters{1}.Distance_InterCellSpacing_Mean_range20]';
CellDistMean = CellDistMean(SeriesID(:,1));
NearestDist = [cellparameters{1}.Distance_ToNearestNeighbor]';
NearestDist = NearestDist(SeriesID(:,1));
NematicOrderParameter = [cellparameters{1}.Architecture_NematicOrderParameter_range30]';
NematicOrderParameter = NematicOrderParameter(SeriesID(:,1));

Alignment_Flow = [cellparameters{1}.Alignment_Flow]';
Alignment_Flow = Alignment_Flow(SeriesID(:,1));
Alignment_Radial = [cellparameters{1}.Alignment_Radial]';
Alignment_Radial = Alignment_Radial(SeriesID(:,1));
Alignment_Zaxis = [cellparameters{1}.Alignment_Zaxis]';
Alignment_Zaxis = Alignment_Zaxis(SeriesID(:,1));

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

%% Difference of parameters
d_dist_surface = dist_surface_2 - dist_surface;
d_LocalDensity = LocalDensity_2 - LocalDensity;
d_LocalNumDensity = LocalNumDensity_2 - LocalNumDensity;
d_CellDistMin = CellDistMin_2 - CellDistMin;
d_CellDistMean = CellDistMean_2 - CellDistMean;
d_NearestDist = NearestDist_2 - NearestDist;
d_NematicOrderParameter = NematicOrderParameter_2 - NematicOrderParameter;

d_Alignment_Flow = Alignment_Flow_2 - Alignment_Flow;
d_Alignment_Radial = Alignment_Radial_2 - Alignment_Radial;
d_Alignment_Zaxis = Alignment_Zaxis_2 - Alignment_Zaxis;

%% Cell position
x = []; y = []; z = [];
for m = 1: length(SeriesID)
    x(m, :) = [cellpos_1(SeriesID(m,1), 2), cellpos_2(SeriesID(m,2), 2)];
    y(m, :) = [cellpos_1(SeriesID(m,1), 1), cellpos_2(SeriesID(m,2), 1)];
    z(m, :) = [cellpos_1(SeriesID(m,1), 3), cellpos_2(SeriesID(m,2), 3)];
end
%% Extract lower positions in biofilm
n_bottom10per = round(length(SeriesID)/20);

[~, idx] = mink(z(:,1), n_bottom10per);
devbottom_x = mean(x(idx,2) - x(idx,1));
devbottom_y = mean(y(idx,2) - y(idx,1));
devbottom_z = mean(z(idx,2) - z(idx,1));

x(:, 2) = x(:, 2) - devbottom_x;
y(:, 2) = y(:, 2) - devbottom_y;
z(:, 2) = z(:, 2) - devbottom_z;

cellpos_2(:,2) = cellpos_2(:,2) - devbottom_x;
cellpos_2(:,1) = cellpos_2(:,1) - devbottom_y;
cellpos_2(:,3) = cellpos_2(:,3) - devbottom_z;


%% Calculate displacements
dx = diff(x,1,2);
dy = diff(y,1,2);
dz = diff(z,1,2);

%% Positions
x1 = x(:,1)-mean(cellpos_1(:,2));
y1 = y(:,1)-mean(cellpos_1(:,1));
z1 = z(:,1)- min(z(:,1));

%% Merge all params in a table
T = table(x1, y1, z1, dist_surface, LocalDensity, LocalNumDensity, NematicOrderParameter...
    , Alignment_Flow, Alignment_Radial, Alignment_Zaxis...
    , dx, dy, dz, d_dist_surface, d_LocalDensity, d_LocalNumDensity, d_CellDistMin, d_CellDistMean, d_NearestDist, d_NematicOrderParameter...
    , d_Alignment_Flow, d_Alignment_Radial, d_Alignment_Zaxis);


end


