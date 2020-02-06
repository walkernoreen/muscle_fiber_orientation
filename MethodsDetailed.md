# Methods: Overview over fiber orientation analysis workflow

## Part 1: Image processing and extraction of fiber orientation measurements
* **Fiji script: fiber_orientation_analysis.py**

### General
* image processing done with a custom jython script written for Fiji
* raw input: A dataset with tiles, zstacks, single channel, czi or tif format
* final output: orientation distribution of muscle fibers in different layers. Format: csv table and visualization images
* other outputs: stitched image, max projected image, flattened image slices along muscle layers
* scripts has GUI where user can define processing settings and steps


### The steps
#### 1. Preprocessing
* Bleaching correction
	* this step: based on a imagej macro by Sean McKinney.
* Downscale the tiles in x & y dimension by factor 0.25
	* purpose: reduce data size for downstream algorithms (surface extraction). Keeping full x&y resolution was not necessary since it did not yield extra information for orientation analysis.
	* important: z-resolution (separation between different muscle layers) must be kept at full resolution.
* Conversion to .tif format

#### 2. Stitching
* images are stitched with Grid/Collection ImageJ plugin: 
	* plugin website: https://imagej.net/Grid/Collection_Stitching_Plugin 
	* publication reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2682522/
* tile ordering: row by row, Right&Down
* also creates a maximum projection for visualization


#### 3. Surface extraction
* surface of planarian is curved, so extracting a single z-slice will not correspond to a single muscle layer but a combination of circular+diagonal+longitudinal layers
* therefore: find the surface of the planarian (i.e. get a height profile of the surface), and extract layers along this surface
* multiple adjacent layers are being extracted: extracted layers 'below' the surface: circular fibers, layers 'above' surface: longitudinal fibers
* extracted layers are then projected to 2D.
* algorithm: ImageJ plugin: Minimum Cost Z Surface Projection:
	* plugin website: MinCostZSurface: https://imagej.net/Minimum_Cost_Z_surface_Projection
	* Important: works for 'rather flat' region in worm center, but not at steep edges (reason: smoothness constraints for surface, and projection to 2D)
* how the plugin works (roughly):
	* algorithm detects the surface along the maximum intensity of the image (implementation actually: finds the minimum intensity ('minimize cost'), therefore invert the image first to create a 'cost image')
	* it has smoothness constraints: surface height cannot change randomly between neighboring pixels (but only by a few z-slices)
* before calling the plugin: image is split into patches (which are larger than the original tiles), for which the surface extraction plugin is called separately.
	* strongly reduces the memory requirements.
	* after surface extraction images are merged again.

#### 4. Derotation/Alignment
* align the worm such that the long axis is aligned vertically.
* done interactively.
* a non-interactive mode is available for development purposes.

#### 5. Orientation analysis
* Orientation is calculated on basis of image intensities
	* Note: Fibers are too dense and too overlapping to segment, therefore no segmentation-based method used.
* Intuition: 
	* To compute orientation at some pixel location: The image gradient in the local neighborhood gives information on orientation. 
	* Gradient at each pixel is always oriented along the steepest decent, for example transition from white fiber to black background
	* The direction vertical to this gradient is the orientation in the local neighborhood, for example pointing along the white fiber.
* ImageJ plugin: OrientationJ http://bigwww.epfl.ch/demo/orientation/ 
* Algorithm (roughly):  
	* uses the gradient structure tensor of the local neighborhood to compute the orientation (instead of gradient since more robust to noise).
	* returns the orientation angles at every image location.
	* Additionally computes the coherency: a measure of how consistently the local neighborhood is oriented (note: this does not dependent on image brightness, which is good since staining efficiency varies). High coherency: strongly orientated local neighborhood. Low coherency: random orientations=isotropic
*  Computes a weighted histogram of orientations for each extracted layer:
	* contribution fo each pixel is weighted by the coherency.
	* histogram is normalized by number of used pixels.
	* user can draw ROI and only this region of interest is taken into account
* additionally visualizes orientation:
	* 1. vector field : orientation is evaluated per image patch and dispayed as (nematic) vector
	* 2. orientation displayed as color
* weighted histograms are saved for further analysis.


## Part 2: Data analysis: analysis of fiber orientation distributions
* **python jupyter notebook: histogram_data_analysis.ipynb.**
	* helper functions: helpers_histogram_analysis.py

### General
* analysis done in python, within a jupyter notebook

### The steps:
#### 1. Analysis for single dataset:
* plotting of orientation histogram for multiple layers (containing different muscle fiber layers) simultaneously
* additionally produces a plot for which each histogram (i.e. per layer) is normalized independently to its max value.
	* allows to visualize diagonal fiber layer, which shows weak coherency and therefore low peak values otherwise.

#### 2. Combined analysis of multiple datasets
* using multiple replicates to obtain statistics
* for each replicate: the layer with the strongest circular resp. longitudinal orientation is detected.
* Then: average these orientation distributions and compute standard deviation.

