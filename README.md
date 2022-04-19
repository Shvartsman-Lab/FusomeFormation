# FusomeFormation

## Requirements

1. MATLAB (versions 2019b, 2021a, and 2022a were used in creating codes)
2. [ilastik](https://www.ilastik.org/) (versions 1.3.0 and 1.3.2 were used for training)
3. FIJI (for image pre-processing)
4. the following MATLAB open-source programs: [munkres.m](https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm), some version of [imagesc3D.m](https://www.mathworks.com/matlabcentral/fileexchange/66638-imagesc3d)

## Getting Started

Clone the repository:
    
    git clone https://github.com/Shvartsman-Lab/FusomeFormation.git
    
## Running Codes

uiopen is a slightly edited version of the MATLAB file that is needed for direct reading of the .h5 files output from ilastik

spec_struct is the main file for reading, reconstructing, assigning proper identity, and measuring the volumes of all fusome pieces.  This script runs split_fusome, allowing for the split channels containing information about fusome and ring canal positional information to be used to split the fusome into its component pieces.  The script then takes these objects and reassigns their labels to match the 2-, 4-, 8-, or 16-cell clusters and then measures the volumes of each piece for further analysis.  The final portion of this script makes a 3D reconstruction of the input data with color labels for each fusome piece, matching the colors of each cell in the network as in previous work (Imran Alsous et al., Current Biology, 2017).

split_fusome uses the input data of the fusome and ring system to recursively split the fusome into its component pieces.  In short, the function cuts the fusome piece at each loop at the location of a chosen ring.  These pieces are then reformed into 2 separate parts, the ring is removed from the selection pool, and then the function repeats itself on one of the two split fusome portions.  This repeats until there are no rings remaining for a given piece of fusome.  From here, the function jumps to another branch of the recursive loop until no rings remain.  The output of this function for n rings are n+1 objects designating the n+1 fusome pieces.

fusome_formation_model compares the measured data of all cluster sizes and compares it to the various volume fractions derived for a range of values of alpha (the volume fraction of fusome in cell 1 after the first division) and beta (the fraction of new fusome material in each plug that remains in the parent cell).  The results of this script create the plots seen in Figure 4.

spec_plots tabulates the data collected in the accompanying Excel file.
