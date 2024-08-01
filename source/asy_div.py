import numpy as np
import pandas as pd

def asy_div(data: np.ndarray, label: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate cytosolic mean GFP intensities for cells with and without foci.

    Args:
        data (np.ndarray): A 2D array where each row represents a cell with columns containing
                           cell area, mean GFP intensity, and foci-related data.
                           - Column 0: Cell area
                           - Column 1: Cell mean GFP intensity
                           - Column 5: Foci area
                           - Column 6: Foci mean GFP intensity
        label (np.ndarray): A 1D array labeling each cell, where each label corresponds to a cell's identifier.

    Returns:
        tuple[np.ndarray, np.ndarray]: Two arrays containing cytosolic mean GFP intensities:
            - GFP0: For cells with no foci.
            - GFP1: For cells with at least one foci.

    Example:
        data = np.array([[1.0, 2.0, np.nan, np.nan, np.nan, 0.5, 3.0],
                         [1.5, 2.5, np.nan, np.nan, np.nan, 0.7, 3.5],
                         [2.0, 3.0, np.nan, np.nan, np.nan, 0.8, 3.8]])
        label = np.array([1, 2, 1])

        GFP0, GFP1 = asy_div(data, label)
        print("GFP0:", GFP0)
        print("GFP1:", GFP1)
    """
    # Find unique labels and corresponding cell indices
    label_unq, i_cells = np.unique(label, return_index=True)
    
    # Extract cell mean GFP and area
    OI23_cell_meanGFP = data[i_cells, 1]
    OI23_cell_area = data[i_cells, 0]
    
    # Calculate total GFP for cells
    OI23_cell_totalGFP = OI23_cell_meanGFP * OI23_cell_area
    
    # Extract foci mean GFP and area, replace NaNs with 0
    OI23_foci_meanGFP = np.nan_to_num(data[:, 6])
    OI23_foci_area = np.nan_to_num(data[:, 5])
    
    # Calculate total GFP for foci
    OI23_foci_totalGFP = OI23_foci_meanGFP * OI23_foci_area
    
    # Initialize arrays for cytosolic mean GFP and total area
    OI23_cytosol_meanGFP = np.zeros(len(i_cells))
    total_area = np.zeros(len(i_cells))
    
    # Loop over each unique cell
    for i, cell_label in enumerate(label_unq):
        # Find all foci indices associated with the current cell
        ind_temp = np.where(label == cell_label)[0]
        
        # Calculate total GFP and area for foci in the current cell
        total_GFP_temp = np.sum(OI23_foci_totalGFP[ind_temp])
        total_area_temp = np.sum(OI23_foci_area[ind_temp])
        
        # Calculate cytosolic mean GFP for the current cell
        if (OI23_cell_area[i] - total_area_temp) != 0:
            OI23_cytosol_meanGFP[i] = (OI23_cell_totalGFP[i] - total_GFP_temp) / (OI23_cell_area[i] - total_area_temp)
        else:
            OI23_cytosol_meanGFP[i] = np.nan  # Handle divide by zero
        
        # Store total area
        total_area[i] = total_area_temp
    
    # Indices where total area of foci is 0 and not 0
    ind0 = total_area == 0
    ind1 = total_area != 0
    
    # Extract GFP0 and GFP1 based on the indices
    GFP0 = OI23_cytosol_meanGFP[ind0]
    GFP1 = OI23_cytosol_meanGFP[ind1]
    
    # Remove NaN values
    GFP0 = GFP0[~np.isnan(GFP0)]
    GFP1 = GFP1[~np.isnan(GFP1)]
    
    return GFP0, GFP1
