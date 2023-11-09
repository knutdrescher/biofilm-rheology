function T = Calculate_Parameters_SingelTrajectories_f1f3(path, param_track_f12, param_track_f23)

path_parametersmat = [path, '\parameters.mat'];
load(path_parametersmat);

path_1to3tracks = [path, '\trackedIDs\all_info.mat'];
index_Series = 3;

load(path_1to3tracks);

%% Extract single cell parameters
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

dist_surface_3 = [cellparameters{3}.Distance_ToSurface_resolution20]';
dist_surface_3 = dist_surface_3(SeriesID(:,index_Series));
LocalDensity_3 = [cellparameters{3}.Architecture_LocalDensity_range30]';
LocalDensity_3 = LocalDensity_3(SeriesID(:,index_Series));
LocalNumDensity_3 = [cellparameters{3}.Architecture_LocalNumberDensity_range30]';
LocalNumDensity_3 = LocalNumDensity_3(SeriesID(:,index_Series));
CellDistMin_3 = [cellparameters{3}.Distance_InterCellSpacing_Min_range20]';
CellDistMin_3 = CellDistMin_3(SeriesID(:,index_Series));
CellDistMean_3 = [cellparameters{3}.Distance_InterCellSpacing_Mean_range20]';
CellDistMean_3 = CellDistMean_3(SeriesID(:,index_Series));
NearestDist_3 = [cellparameters{3}.Distance_ToNearestNeighbor]';
NearestDist_3 = NearestDist_3(SeriesID(:,index_Series));
NematicOrderParameter_3 = [cellparameters{3}.Architecture_NematicOrderParameter_range30]';
NematicOrderParameter_3 = NematicOrderParameter_3(SeriesID(:,index_Series));

Alignment_Flow_3 = [cellparameters{3}.Alignment_Flow]';
Alignment_Flow_3 = Alignment_Flow_3(SeriesID(:,index_Series));
Alignment_Radial_3 = [cellparameters{3}.Alignment_Radial]';
Alignment_Radial_3 = Alignment_Radial_3(SeriesID(:,index_Series));
Alignment_Zaxis_3 = [cellparameters{3}.Alignment_Zaxis]';
Alignment_Zaxis_3 = Alignment_Zaxis_3(SeriesID(:,index_Series));

%% Difference of parameters
d_dist_surface = dist_surface_3 - dist_surface;
d_LocalDensity = LocalDensity_3 - LocalDensity;
d_LocalNumDensity = LocalNumDensity_3 - LocalNumDensity;
d_CellDistMin = CellDistMin_3 - CellDistMin;
d_CellDistMean = CellDistMean_3 - CellDistMean;
d_NearestDist = NearestDist_3 - NearestDist;
d_NematicOrderParameter = NematicOrderParameter_3 - NematicOrderParameter;

d_Alignment_Flow = Alignment_Flow_3 - Alignment_Flow;
d_Alignment_Radial = Alignment_Radial_3 - Alignment_Radial;
d_Alignment_Zaxis = Alignment_Zaxis_3 - Alignment_Zaxis;

%% Calculate displacements
dx = param_track_f12.dx + param_track_f23.dx;
dy = param_track_f12.dy + param_track_f23.dy;
dz = param_track_f12.dz + param_track_f23.dz;

%% Positions
x1 = param_track_f12.x1;
y1 = param_track_f12.y1;
z1 = param_track_f12.z1;

%% Merge all params in a table
T = table(x1, y1, z1, dist_surface, LocalDensity,LocalNumDensity, NematicOrderParameter...
    , Alignment_Flow, Alignment_Radial, Alignment_Zaxis...
    , dx, dy, dz, d_dist_surface, d_LocalDensity, d_LocalNumDensity, d_CellDistMin, d_CellDistMean, d_NearestDist, d_NematicOrderParameter...
    , d_Alignment_Flow, d_Alignment_Radial, d_Alignment_Zaxis);

end


