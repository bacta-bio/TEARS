# Segementations of E. coli with synthetic organelles

Supporting code for research article: Spatial engineering of E. coli with addressable phase-separated RNAs

Cellpose, a generalist algorithm for cellular segmentation [Stringer C et al. 2021] is trained and then segments E. coli cells from 16-bit grayscale microscopy images, resulting in binary mask images. Segmentation errors are corrected semi-automatically by an in-house ImageJ script (shared code,  segmentation_correction.zip). The corrected segmentations are then re-used as input to the training dataset to increase the accuracy of the training model. A coarse split or merge by the painting tools by ImageJ is used to correct the segmentation and the final segmentation is then achieved by re-running the script. The Regions of Interest (ROI) of each cell are extracted and the relevant information as intensity, area are output into a csv file.


For the droplet fusion analysis, droplets were first segmented by our ImageJ script (shared code, fusion_display.zip) to generate ROIs (Region of Interest) which are saved into a zip file. By using the script, the relevant information of each ROI (Circularity, Area, Intensity) of each frame is extracted and exploited into a csv file, the csv file is then read and interpreted by the shared python script (shared code, droplets_info_extraction.zip) to generate the fusion process display. Examples are included in the shared code. 
