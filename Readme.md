
## Summary

This project analyzes the orientation of muscle fibers in planarian flatworms. 

The scripts were developed for the Rink Lab at MPI-CBG in Dresden.

Scripts are written by Noreen Walker, Scientific computing facility at MPI-CBG. </br>
The bleaching correction step is based on a Macro by Sean McKinney.

## Usage
1) First run `fiber_orientation_analysis.py` in Fiji. 
2) Then run the jupyter notebook `histogram_data_analysis.ipynb.`

### Usage details step 1:
* Open the script in Fiji, and click on *Run* in the script editor.
* A GUI with these **parameters** will show up:
	- *Input folder*: The file containing your raw tiles (with/without masterfile).
	- *Results folder*: Root directory for results. A subdirectory with the name of the dataset is then automatically created.
	- *File type raw data*: czi or tif
	- *Do background (bleaching) correction*: correct for inhomogeneous brightness caused by bleaching at the tile-overlap region
	- *Downscale raw tile in x and y direction*: keep checked when using the full pipeline.
	- *Number of tiles in x resp. y direction*: input for stitching
	- *Number of slices to be extracted above/below fitted surface*: Try ca 10.
	- *Derotation mode*: choose between no/automatic/user-defined horizontal alignment
	- *ROI determination mode*: choose between no/automatic/user-defined determination of region in which the orientation histogram is calculated
	- *Grid size for vectorfield patches*: for final vectorfield orientation visualization
	- *Surface extraction: image patch size*: splits image into patches of this size before extracting the layers. saves memory
	- *Processing steps*: all, first part (automatic), second part (interactive), custom
	- *First/Last step to run*: only if 'custom' mode: specify which steps to run

### Usage details on step 2:
* Usage is described in the notebook itself.
* For an example output see *histogram_data_analysis_for_illustration.html*


## Example data
* A **small example dataset** which contains a subset of tiles (and lower z-sampling) can be found in [here](https://github.com/walkernoreen/muscle_fiber_orientation).
* Run the Fiji script with the following parameters: 
	* *File type raw data*: tif
	* *Do background correction*: checked
	* *Downscale raw tiles*: checked
	* *Number of tiles in x resp. y direction*: 3 resp. 1
	* *Surface extraction: number of slices above/below*: 5
	* *Derotation mode*: 2
	* *ROI determination mode*: 2
	* *Grid size*: 50 (default)
	* *Surface extraction: image patch size*: 500 (default)
	* *Processing steps*: all

## Software requirements:
### Fiji
* Activate the following update site (*Help -> Update... -> ManageUpdateSites*):
	* `OrientationJ`
* Download the `MinCostZSurface` jar from https://imagej.net/Minimum_Cost_Z_surface_Projection, then install via *Plugins -> InstallPlugin*
### Python
* Any standard python3 installation will do.

## Required input data format
* Raw input is a collection of z-stack tiles (czi or tif) displaying fluorescently labelled muscle fibers, single channel.


## Workflow description
For a more detailed explanation of the methods please check the [methods](MethodsDetailed.md) document.

### Part 1: fiber_orientation_analysis.py
* **From raw data to orientation measurements.**
* Input: Raw input is a collection of z-stack tiles (czi or tif) displaying fluorescently labelled muscle fibers, single channel.
* Output: 
	* data: csv tables with orientation measurements for each sliced muscle layer, additionally plotted as histograms.
	* images: stitched z-stack, zmax projection, image of flattened-surface (muscle layers), muscle layers colored by fiber orientation & vectorfield of orientation
* **Main processing steps** are:
	1. downscale individual tiles in x&yand correct for bleaching
	2. stitch tiles (Grid/Collection Stitching)
	3. extract and flatten surface (MinCostZSurface).
	4. derotate the stitched image (align vertically)
	5. do intensity based orientation analysis (OrientationJ). 

### Part 2: histogram_data_analysis.ipynb
* **analyze orientation histograms**
* plots orientation distribution for each layer, and can combine results of multiple experiments.

## Tips
* Memory requirements for extract_surface & orientation_analysis are rather high when processing a full dataset.
	* Use the option to split image into patches for surface extraction.
	* Use a workstation if needed
	* Make sure that maximum allowed memory of Fiji is set to high value.
* Do not use blank spaces in file names.
* Downstream analysis steps require specific filenames, produced by upstream steps. Changing the filenames of intermediate steps will likely break the workflow.


