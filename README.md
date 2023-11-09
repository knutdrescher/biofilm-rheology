# Cell tracking for continuous cell displacements based on single-cell segmented labels

Code was written in Matlab R2020b

The code provides 3D cell tracking.
This algorithm is under assumption as follows:

- Cells are continuously displaced 
- Neighbor cells should keep locational relationships
- Cell locations should not be switched locally between time-frames

## Version 1.0 (28.09.2023)

**This version is optimized only for dataset in our publication (*Ohmura, et al., submitted*, 2023).**\
The dataset describes deformation of mesoscaled bacterial biofilms under shear flow. The tracked dataset is available in our Zenodo repository (URL).

## Input data

- Time series of 3D voxels of cell labels

To obtain cell labels, apply cell segmentation for microscope cell images.
In case of bacterial biofilms, we recommend using StarDist OPP, which is specialized for single-cell segmentation of biofilms (Jelli, Ohmura, Netter, et al., *Mol. Microbiol.*, **119**, 659, 2023)

## Usage of the scripts

### 1. Extract single-cell parameters

Extract single-cell parameters from images with cell labels by using BiofilmQ (Hartmann, et al., *Nat. Microbiol.*, **6**, 151, 2021).

- Download BiofilmQ (<https://drescherlab.org/data/biofilmQ/docs/>), which is a GUI for cell segmentation and calculating single-cell parameters
- See a section "Segmentation > Import segmentation" (Label image input)
- After importing cell labels, see a section "Calculate parameters"

article (*Ohmura, et al., submitted*, 2023), this process corresponds to Step 1 in Figure S1.

### 2. Track a cell manually

Track a cell manually as a reference. Put the tracked label into a MAT file.\
In the article (*Ohmura, et al., submitted*, 2023), this process corresponds to Step 2 in Figure S1.

### 3. Track all cells automatically

Execute "TrackInAllFrames.m" with the local paths to the dataset produced by BiofilmQ and the manual-tracked label.\
In the article (*Ohmura, et al., submitted*, 2023), this process corresponds to Step 3 in Figure S1.

### 4. Calculate cell-tracked trajectories and orientation differences

Execute "CalculateParameters_From_TrackingIDs.m" with the local path to the dataset produced by BiofilmQ
