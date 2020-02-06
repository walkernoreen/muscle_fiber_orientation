# ==================================================================================
# FOR DETAILS AND INSTRUCTIONS ON USAGE: see Readme.md
# ==================================================================================
#
# This script analyzes the orientation of muscle fibers in planarian flatworms.
# It does the following steps: preprocessing of individual tiles, stitching, surface extraction, vertical alignment, intensity based orientation analysis.
#
# Input: a folder containing multiple tiles of single channel z-stacks of planarians. format is .czi or .tif.
# Output: Visualization of fiber orientation images, .csv files with orientation information. The .csv files serve as input for downstream analysis with histogram_data_analysis.ipynb
#
# Author: Noreen Walker, Scientific Computing Facility, MPI-CBG
# 		  The bleaching correction step is based on a Macro by Sean McKinney.
#
# 2020-02-10: v1.4.2 first public version.
#				 updated code documentation and removed obsolete functions. otherwise identical to internal version 1.4.1.
#
# ==================================================================================
# Important for developers:
# *	Keep the string choices in processing_steps_type, custom_first_step, custom_last_step, the function selectProcessingSteps() and the substep selection criteria in main() strictly in sync!
#   Also keep them loosely in sync with the output directory names.
# * Part of OrientationJ Vectorfield was reimplemented since spurious error that sometimes only first slice was plotted
# ==================================================================================


#@ String (visibility=MESSAGE, value="version 1.4.2") msg
#@File (label= "Input folder (raw data)" , style= "directory") inputdir_java
#@File (label= "Results folder (root folder for all results)" , style= "directory") resultsrootdir_java
#@String (label="File type raw data", choices={"czi", "tif"}) extension
#@Boolean (label="Do background (bleaching) correction",value=True) do_background_correction
#@Boolean (label="Downscale raw tiles in x and y direction by factor 0.25 (for full pipeline: check box, for stitching only: uncheck box)",description="z direction is never downscaled. note: full image must fit into RAM memory. The full pipeline can probably also run without downscaling but memory requirements will be large and benefit is doubtful",value=True) do_xy_downscaling
#@Integer (label="Number of tiles in x direction") numtiles_x
#@Integer (label="Number of tiles in y direction") numtiles_y
#@Integer (label="Surface extraction: Number of slices extracted above the fitted surface (circular fibers)", description="Total number of extracted slices: num_slices_above+numslices_below+1 (centerslice). Note: final value can deviate when using surface extraction in patches.",value=11) numslices_above
#@Integer (label="Surface extraction: Number of slices extracted below the fitted surface (longitudinal fibers)", description="Total number of extracted slices: num_slices_above+numslices_below+1 (centerslice). Note: final value can deviate when using surface extraction in patches.",value=11) numslices_below
#@Integer (label="Derotation mode for vertical alignment (0: only rough alignment (90deg), 1: automatic (very basic), 2: interactive by user)",value=1) derotation_mode
#@Integer (label="ROI determination mode for orientation histogram (0: use full image, 1: automatic (very basic), 2: interactive by user)",description="determines region over which orientation histogram is calculated", value=1) mask_mode
#@Integer (label="Grid size of patches for orientation vectorfield visualization (default: 50 px)", description="patches are square with length=width=this value", value=50) vector_grid_size
#@Integer (label="Surface extraction: Size of patches (in pixel) for surface extraction (net size). default=500",description="Stitched image is split into patches before layer extraction to save memory. The actual patch_size is ca 20% larger because an overlap with neighboring patches is added. Lower value for patch size: Faster and less memory required, but lower accuracy. If set to -1: use full image (single patch)",value=500) patch_increment
#@String (label="Processing steps", choices={"all", "first part (all automatic steps: up to stitching (incl))", "second part (all interactive steps: from derotation on)","custom"}) processing_steps_type
#@String (label="Only if 'custom' mode for processing steps: First step to run",description="preprocess: downscale & correct for uneven background", choices={"preprocess","stitch","extract surface","derotate","orientation analysis"},style="radioButtonHorizontal") custom_first_step
#@String (label="Only if 'custom' mode for processing steps: Last step to run",description="preprocess: downscale & correct for uneven background", choices={"preprocess","stitch","extract surface","derotate","orientation analysis"},style="radioButtonHorizontal") custom_last_step
#@OpService op



from ij import IJ
from loci.plugins import BF
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator
from net.imglib2.img.display.imagej import ImageJFunctions
from ij.plugin import FolderOpener
from ij.plugin import ZProjector
from loci.formats import ImageReader
from loci.plugins.in import ImporterOptions
from ij.measure import Measurements
from ij.gui import WaitForUserDialog
from ij import WindowManager
from ij.plugin import Macro_Runner
from ij.gui import Overlay
from ij.gui import Plot
from ij.plugin import ImageCalculator
from ij.plugin import HyperStackConverter
import os
import sys
import glob
import collections
import math
import time
from ij.measure import ResultsTable
from ij.gui import Roi, PolygonRoi, Line
from java.awt import Color
from string import digits




# ==================================================================================
# ======= general parameters + paths =========
#GUI: rawdata_dir, results_rootdir, extension
downscaled_subdir="1_preprocessing"
stitchprep_subdir="2a_tmp_stitch"
stitchfinal_subdir="2b_stitched"
surfaceextract_subdir="3_surface_extraction"
derotated_subdir="4_derotated"
orientation_subdir="5_orientation_analysis"


# ===== params step 1: background correction & downscale in xy ======
if do_xy_downscaling: # see GUI
	scalefactor_xy=0.25 # between ]0,1]
else:
	scalefactor_xy=1


# ===== params step 2: stitching ======
#GUI: numtiles_x, numtiles_y
title_fused_image="full_zstack_stitched" # without the tif

if do_xy_downscaling:
	min_r=0.7 # 0.7 default # pairwise correlation filtering (x,y,z options: for now hardcoded below) (keep cst)
else:
	min_r=0.5
tiles_overlap=15 # in (%), for x and y as initial alignment (keep cst)


# ==== params step 3: surface extraction of tiles ====
downsamplefactor_xy=0.25 # mincostzsurface: downsampling for surface extraction.
downsamplefactor_z=1 # mincostzsurface: downsampling for surface extraction.
max_dz=5 # mincostzsurface: max change in z for neighboring pixels (creates smoothness constraint)
sigma_cost=0 # (Deprecated).gaussian blur of the cost function (zstack) to extract a smoother surface (avoid muscle fiber pattern). (active if >0).
sigma_zmap=11 # gaussian blur of the calculated height map to extract a smoother surface (avoid muscle fiber pattern). (active if >0) # NOTE: Applying 2 sigmas could be redundant ...
percentage_patch_overlap=15 # resave the image in patches: how much overlap. The true patch size is patch_increment/(1-percentage_patch_overlap/100)

# ==== params step 4: derotation ====
#GUI: derotation_mode

# ==== params step 5: orientation analysis
flattenedsurface_filename_contains= "flattened_surface" # exclude the full stack
flattenedsurface_filename_contains_not="zmax"
# GUI: mask_mode


# ==================================================================================
# main processing function 
# ==================================================================================
def main():

	# clean up from previous round
	IJ.run(None, "Close All", "") 

	# extract which steps to process
	steps=selectProcessingSteps(processing_steps_type, custom_first_step, custom_last_step)
	IJ.log("\n *** Started processing *** ")
	IJ.log("\nWill execute processing steps: "+str(steps)+"\n")

	# convert java paths to strings
	rawdata_dir=inputdir_java.getPath()
	results_rootdir=resultsrootdir_java.getPath()

	# infer the dataset name:
	fs=os.sep
	dataset_name=rawdata_dir.split(fs)[-1] 
	IJ.log("Processing dataset with (inferred) name: "+dataset_name)

	# create names of save directories
	parent_dir=os.path.join(results_rootdir,dataset_name) # shared parentdir
	downscaled_dir=os.path.join(parent_dir,downscaled_subdir)
	surfaceextract_dir=os.path.join(parent_dir,surfaceextract_subdir)
	stitchprep_dir=os.path.join(parent_dir,stitchprep_subdir)
	stitchfinal_dir=os.path.join(parent_dir,stitchfinal_subdir)
	derotated_dir=os.path.join(parent_dir,derotated_subdir)
	orientation_dir=os.path.join(parent_dir,orientation_subdir)
	
	
	# step 1: background (bleaching) correction and downscaling
	if "preprocess" in steps:
		tstart=time.time()
		if do_background_correction:
			background_imp=computeBackgroundImage(rawdata_dir,downscaled_dir, extension) # background_imp is full orig size
		else:
			background_imp=None
		preprocessImages(rawdata_dir,downscaled_dir, extension, scalefactor_xy, background_imp)
		print "Time: preprocessing took        "+ str(int(time.time()-tstart))+ " sec."	


	# step 2: stitching
	if "stitch" in steps:
		tstart=time.time()
		stitchImagesGridCollection(downscaled_dir,stitchprep_dir, stitchfinal_dir, numtiles_x, numtiles_y,
			min_r,tiles_overlap, extension="tif", savename=title_fused_image,do_zmax_projection=True)
		print "Time: stitching took            "+ str(int(time.time()-tstart))+ " sec."


	# step 3: extract flattened surfaces of tiles
	if "extract surface" in steps:
		tstart=time.time()
		extractFlattenedSurfaces(stitchfinal_dir, surfaceextract_dir, downsamplefactor_xy, downsamplefactor_z,
								 max_dz, numslices_above, numslices_below, sigma_cost, sigma_zmap, extension="tif", patch_increment=patch_increment)
		print "Time: surface extraction took   "+ str(int(time.time()-tstart))+ " sec."


	# step 4: derotation (vertical alignment: head at top)
	if "derotate" in steps:
		tstart=time.time()
		derotateImages(stitchfinal_dir, surfaceextract_dir, derotated_dir, derotation_mode, zmax_identifier="zmax",flatsurface_identifier=flattenedsurface_filename_contains, extension="tif")
		print "Time: derotation took           " + str(int(time.time()-tstart))+ " sec. (incl optional wait for user action)"	
		

	# step 5: orientation analysis
	if "orientation analysis" in steps:
		tstart=time.time()
		analyzeFiberOrientation(derotated_dir, orientation_dir, flattenedsurface_filename_contains,flattenedsurface_filename_contains_not, mask_mode, vector_grid_size, extension="tif")
		print "Time: orientation analysis took "+ str(int(time.time()-tstart))+ " sec. (incl optional wait for user action)"

	IJ.log("\n *** Done! *** ")



# ==================================================================================
# preparation
# ==================================================================================
def selectProcessingSteps(processing_steps_type, custom_first_step, custom_last_step):
	"""Processing the user selection of which steps to execute and returns a list of named steps.
	Important to keep up-to-date with user-input (@Parameters) and selection criteria in main()"""

	all_steps_in_order=["preprocess","stitch","extract surface","derotate","orientation analysis"]
	
	#processing_steps_type can be: "all", "first part (automatic steps)", "second part (interactive steps)","custom"
	
	if processing_steps_type=="all":
		steps=all_steps_in_order
	elif processing_steps_type=="first part (all automatic steps: up to stitching (incl))":
		steps=all_steps_in_order[:3] # preprocess, extract surface, stitch
	
	elif processing_steps_type=="second part (all interactive steps: from derotation on)":
		steps=all_steps_in_order[3:] # derotate, orientation analysis
	
	elif processing_steps_type=="custom":
		idx0=all_steps_in_order.index(custom_first_step)
		idx1=all_steps_in_order.index(custom_last_step)
		if idx1<idx0: # wrong user input
			steps=[]
		else:
			steps=all_steps_in_order[idx0:idx1+1]	
			
	return steps


# ==================================================================================
# background correction and downscaling step: main functionality
# ==================================================================================
def computeBackgroundImage(datadir, savedir, extension):
	"""
	Computes the average background image. Loops through all images in datadir and computes the backgroudn image by
	doing z-rpojection, then strong blurring (sigma ca 10% of image size), then averaging. The resulting image is 
	returned and, for reference, saved to disk.
	Input images can be czi (with/without masterfile), tif, etc. (defined by extension). 
	Quite similar in structure to preprocessImages().
	Returns: background_imp: 2D image of flatfield, normalized to mean=1, original image widhth x height. 
							Use it for background/bleaching correction via division.
	"""
	IJ.log("\nStarting computeBackgroundImage ...")

	verbose=True

	assert datadir!=savedir
	assert os.path.exists(datadir)

	# initialize save dir (actually: a new subdir)
	savedir=os.path.join(savedir,"background_image");
	if not os.path.exists(savedir):
		os.makedirs(savedir)

	# check if masterfile exists
	folder_has_masterfile, masterfile=searchForMasterfile(datadir,extension,False)

	# initialize collection of background images of all tiles
	bg_impstack=None

	# image loading depends on yes/no existence of masterfile
	if folder_has_masterfile:
		# check how many files belong to the series
		bfreader=ImageReader()
		bfreader.setId(masterfile)
		seriescount=bfreader.getSeriesCount()
		if verbose:
			IJ.log("Nr of images belonging to masterfile (seriescount): "+ str(seriescount))

		# loop over images
		for idx in range(1,seriescount+1): # one-based
			if verbose:
				tstart=time.time()

			# load image (can be significantly slower than without masterfile)
			IJ.log("Processing file with idx: "+ str(idx))
			imp=loadImagesViaMasterfile(masterfile,idx,verbose=False)[0]

			# create background image
			bg_imp=backgroundImageFromSingleTile(imp)
	
			# collect all background images in a stack
			if bg_impstack is None: # first iteration
				bg_impstack=bg_imp	
			else:
				stack = bg_impstack.getStack() 
				stack.addSlice(bg_imp.getProcessor())		
				bg_impstack.setStack(stack); 

	# no master file exists
	else:
		# get list with all imagefiles in data folder
		imagefiles=getFileList(datadir,extension=extension, verbose=verbose)
		
		# intialize import options
		bfoptions=ImporterOptions()
		
		# loop over all images
		for idx in range(len(imagefiles)): # zero-based
			if verbose:
				tstart=time.time()

			bfoptions.setId(imagefiles[idx])
			IJ.log("Processing: idx="+str(idx)+ ": "+ imagefiles[idx])
		
			# load image
			imp=BF.openImagePlus(bfoptions)[0]

			# create background image
			bg_imp=backgroundImageFromSingleTile(imp)
	
			# collect all background images in a stack
			if bg_impstack is None: # first iteration
				bg_impstack=bg_imp	
			else:
				stack = bg_impstack.getStack() 
				stack.addSlice(bg_imp.getProcessor())		
				bg_impstack.setStack(stack); 


	# average the background images
	background_imp = ZProjector.run(bg_impstack,"avg")
	meanval=background_imp.getStatistics().mean
	IJ.run(background_imp, "Divide...", "value="+str(meanval));

	# save the images to disk for debugging
	IJ.save(bg_impstack,os.path.join(savedir,"background_stack_unnormalized.tif"))
	IJ.save(background_imp,os.path.join(savedir,"background_image_origsize.tif"))

	IJ.log("Finished computation of background image.\n")

	return background_imp



def preprocessImages(datadir, savedir, extension, scalefactor_xy, background_imp=None):
	"""
	Loops through all images in datadir, downscales them in x and y by scalefactor_xy (if <1), and saves them in savedir.
	Input images can be czi (with/without masterfile), tif, etc. (defined by extension). 
	Corrects for background image (by division) if background_imp is not None.
	"""
	IJ.log("\nStarting preprocessImages ...")
	if scalefactor_xy<1:
		IJ.log(".. .will downscale images in x&y by factor "+str(scalefactor_xy));
	else:
		IJ.log("... will not downscale images.")

	verbose=True

	assert scalefactor_xy>0
	assert datadir!=savedir
	assert os.path.exists(datadir)

	# initialize save dir
	if not os.path.exists(savedir): 
		os.makedirs(savedir)

	# check if masterfile exists
	folder_has_masterfile, masterfile=searchForMasterfile(datadir,extension,False)

	# image loading depends on yes/no existence of masterfile
	if folder_has_masterfile:
		# check how many files belong to the series
		bfreader=ImageReader()
		bfreader.setId(masterfile)
		seriescount=bfreader.getSeriesCount()
		if verbose:
			IJ.log("Nr of images belonging to masterfile (seriescount): "+ str(seriescount))

		# loop over images
		for idx in range(1,seriescount+1): # one-based
			if verbose:
				tstart=time.time()

			# load image (can be significantly slower than without masterfile)
			IJ.log("Processing file with idx: "+ str(idx))
			imp=loadImagesViaMasterfile(masterfile,idx,verbose=False)[0]

			# create savename
			orig_title=imp.getTitle() 
			savename="rescaled_"+orig_title.split(".",1)[0]+"("+str(idx)+").tif"
			
			# process image
			RescalingBackgroundCorrectionAndSaving(imp,scalefactor_xy, savedir, savename, background_imp)

	else:
		# get list with all imagefiles in data folder
		imagefiles=getFileList(datadir,extension=extension, verbose=verbose)
		
		# intialize import options
		bfoptions=ImporterOptions()
		
		# loop over all images
		for idx in range(len(imagefiles)): # zero-based
			if verbose:
				tstart=time.time()

			bfoptions.setId(imagefiles[idx])
			IJ.log("Processing: idx="+str(idx)+ ": "+ imagefiles[idx])
		
			# load image
			imp=BF.openImagePlus(bfoptions)[0]

			# create savename
			orig_title=imp.getTitle()
			savename="rescaled_"+orig_title.rsplit(".",1)[0]+".tif"

			# process image
			RescalingBackgroundCorrectionAndSaving(imp,scalefactor_xy, savedir, savename, background_imp)
			
	IJ.log("Finished downscaling.\n")


# ==================================================================================
# background correction and downscaling step: helper functions
# ==================================================================================
def backgroundImageFromSingleTile(imp):
	""" computes a background image from a single image tile by doing zmax projection, then blurring.
	Background images from all tiles should be averaged in a later step.
	params:
	imp: w x h x z x 1channel x 1frame
	"""
	# z-project
	bg_imp = ZProjector.run(imp,"max");
	IJ.run(bg_imp, "32-bit", "");
	
	# blur image
	# use sigma roughly 10% of the image dimensions (heuristic)
	width,height,_,_,_=bg_imp.getDimensions() 
	sigma=int(0.5*(width+height)/10)
	IJ.run(bg_imp, "Gaussian Blur...", "sigma="+str(sigma))

	return bg_imp


def RescalingBackgroundCorrectionAndSaving(imp, scalefactor_xy, savedir, savename, background_imp=None):
	"""Helper for preprocessImages.
	Downscales stack (imp) in x,y by scalefactor_xy, then saves the stack to disk.
	"""
	# do background correction
	if background_imp is not None:
		IJ.log("Doing background correction.")
		imp = ImageCalculator().run("Divide create 32-bit stack", imp, background_imp);

	# downscale image in x & y
	if scalefactor_xy<1:
		orig_width,orig_height,_,nslices,_=imp.getDimensions() 
		new_width=int(orig_width*scalefactor_xy)
		new_height=int (orig_height*scalefactor_xy)
		IJ.run(imp, "Scale...", "x="+str(scalefactor_xy)+" y="+str(scalefactor_xy)+" z=1.0 width="+str(new_width)+" height="+
			str(new_height)+" depth="+str(nslices)+	" interpolation=Bilinear average process create")
		imp=IJ.getImage() 
	
	# save image as 16 bit (without rescaling, conversion only needed if downscaled)
	if imp.getType()!=ImagePlus.GRAY16:
		IJ.run("Conversions...", " ");
		IJ.run(imp, "16-bit", "");
		IJ.run("Conversions...", "scale"); # reset to default conversion behaviour

	IJ.save(imp,os.path.join(savedir,savename))

	# close windows
	IJ.run(imp, "Close All", "")


				
def searchForMasterfile(datadir,extension,verbose=True):
	"""Helper for preprocessImages.
	Searches directory for a masterfile: defined by the lack of brackets "(" and ")".
	Only czi can have masterfile, if tif then there is never a masterfile.
	"""
	if (extension!="czi" and extension!=".czi"):
		IJ.log("Data is not czi format. Skipped looking for masterfile.")
		return False, None
		
	# get list with all imagefiles in data folder
	imagefiles=getFileList(datadir,extension=extension, verbose=verbose)

	# search for masterfile (i.e. filename without brackets)
	masterfile=None
	for fn in imagefiles:
		basename=fn.split(os.sep)[-1]
		if ("(" not in basename) and (")" not in basename):
			masterfile=fn
			break
	if masterfile is not None:
		IJ.log("Found masterfile: "+ masterfile)
		folder_has_masterfile=True
	else:
		IJ.log("No masterfile found.")
		folder_has_masterfile=False

	return folder_has_masterfile, masterfile
	

def loadImagesViaMasterfile(masterfile, imgidx, verbose=False):
	"""
	Helper for preprocessImages.
	Loads one or more images/stacks from an image folder which is controlled by a masterfile.
	The individual image files are 'series elements' of the master file. Therefore, direct loading of an image would not work
	(would always load the first file of the folder).	
	Note that loading an individual image via a masterfile is slow. Loading multiple image in one go (provide a list) is much more efficient.
	This fct was written for loading files in .czi format (not tested for other formats).
	Params:
	masterfile: full path, e.g. "/somefolder/egfp.czi" (image files are then "egfp(1).czi" etc)
	imgidx: integer or list of ints. one-based(!). which images of series to load
	returns: the list of loaded image(s): output of  BF.openImagePlus(options)	
	"""
	if verbose:
		import time
		tstart=time.time()

	# check how many files belong to the series
	bfreader=ImageReader()
	bfreader.setId(masterfile)
	seriescount=bfreader.getSeriesCount()
	
	# intialize import options
	bfoptions=ImporterOptions()
	bfoptions.setId(masterfile)

	def setToOn(bfoptions,idx, seriescount):
		if idx in range(1,seriescount+1):
			# shift from one to zero based (series ids are here zero-based for some strange reason)
			bfoptions.setSeriesOn(idx-1,True)
		else: 
			IJ.log ("Not a valid image index: "+str(idx)+ ". Skipping ...")
		return bfoptions
		
	# specify the images to be loaded
	if isinstance(imgidx,collections.Iterable):
		for iidx in imgidx:
			bfoptions=setToOn(bfoptions,iidx,seriescount)
	else:
		bfoptions=setToOn(bfoptions,imgidx,seriescount)

	# load image(s)
	imparray=BF.openImagePlus(bfoptions)
	if verbose:
		IJ.log("Elapsed time for loading idx(s) "+ str(imgidx)+": "+str(int(time.time()-tstart))+" sec.")
	return imparray




# ==================================================================================
# stitching step: main functionality
# ==================================================================================
def stitchImagesGridCollection(input_dir, prepstitch_dir,finalsave_dir, numtiles_x, numtiles_y,  min_r=0.7,tiles_overlap=15, extension="tif", savename="fused_image",
							   do_zmax_projection=True):
	"""Stitches all images in input_dir, using Grid/Collection. Intermediate results (slices of stitched iamge) are saved in prepstitch_dir, 
	the TileConfiguration is saved in the input_dir (cannot be changed...), the final stitched image stack in finalsave_dir. 
	Additionally saves a zmax projection of the stitched image.
	Params:
	numtiles_x: tile configuration
	numtiles_y: tile configuration
	min_r: pairwise correlation filtering (x,y,z options: for now hardcoded below)
	tiles_overlap: in (%), for x and y as initial alignment (keep cst)
	extension: "tif" etc
	savename: basename without extension
	do_zmax_projection: whether to additionally save a zmax projection of the stitched image
	"""	
	IJ.log("\nStarting stitchImages (Grid/Collection)...")

	# prepare folders
	assert os.path.exists(input_dir)
	if not os.path.exists(prepstitch_dir): 
		os.makedirs(prepstitch_dir)
	if not os.path.exists(finalsave_dir):
		os.makedirs(finalsave_dir)

	# extract spatial calibration
	IJ.log("Extracting calibration")
	fntmp=getFileList(input_dir, extension)[0]
	tmpImp=BF.openImagePlus(fntmp)[0]
	calib=tmpImp.getCalibration()
	del tmpImp


	# (a) stitching
	IJ.log("Started Grid/Collection Stitching (very slow)")
	
	filename_pattern=getFilenamePatternOfTiles(input_dir,extension) # a string like "tile_({i}).tif"

	# all tiles must fit into memory
	IJ.run(None, "Grid/Collection stitching", "type=[Grid: row-by-row] order=[Right & Down                ] grid_size_x="+str(numtiles_x)+ \
	" grid_size_y="+str(numtiles_y)+" tile_overlap="+str(tiles_overlap)+" first_file_index_i=1 directory=["+input_dir+ \
	"] file_names=["+filename_pattern+"] output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] "+\
	"regression_threshold="+str(min_r)+" max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap ignore_zstage"+\ # +-ignore_zstage
	" computation_parameters=[Save memory (but be slower)] image_output=[Write to disk] output_directory=["+prepstitch_dir+"]")
	# optional: "use_virtual_input_images"

	
	# (b) combine individually saved slices into a stack + resave
	IJ.log("Resaving stitched slices as stack")
	file_identifier="img_"
	imp = FolderOpener.open(prepstitch_dir, " file="+file_identifier+" virtual")

	# recover calibration
	imp.setCalibration(calib);
			

	IJ.saveAs(imp, "Tiff", os.path.join(finalsave_dir,savename+".tif"))
	IJ.log("Finished Grid/Collection Stitching")
	
	# (c) create z-max projection of fused image + save
	if do_zmax_projection:
		IJ.log("doing zmax projection")
		fn=os.path.join(finalsave_dir,savename+".tif")
		IJ.log("...of file"+fn)
		IJ.run("Bio-Formats", "open="+fn+" color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
		imp=IJ.getImage()
		imp = ZProjector.run(imp,"max");
		IJ.save(imp,os.path.join(finalsave_dir,"zmax_of_"+savename+".tif")) 

	IJ.run(imp, "Close All", "")

	IJ.log("Finished stitching (Grid/Collection).\n")



# ==================================================================================
# stitching step: helper functions
# ==================================================================================
def getFilenamePatternOfTiles(input_dir, extension):
	"""creates a string that can used as file_name expression for grid collection plugin. 
	For example tile_({i}).tif if the tiles are called tile_(1).tif, tile_(2).tif.
	Note: assumes that the counter is 1,2,3,... , i.e. not 01,02,03.. or 4,5,6,...
	"""
	fullfiles=getFileList(input_dir,extension)
	files=[] # basenames (incl extension)
	for i in range(len(fullfiles)):
		fn=fullfiles[i].split(os.sep)[-1]
		files.append(fn)

	def longest_common_suffix(list_of_strings):
		# longest commen suffix in a list fo strings
		reversed_strings = [''.join(reversed(s)) for s in list_of_strings]
		reversed_lcs = os.path.commonprefix(reversed_strings)
		lcs = ''.join(reversed(reversed_lcs))
		return lcs

	prefix=os.path.commonprefix(files)
	suffix=longest_common_suffix(files)
	
	filename_pattern=prefix+"{i}"+suffix # grid collection syntax
	return filename_pattern



def getFilenamePatternOfTilesXY(input_dir, extension):
	"""creates a string that can used as file_name expression for grid collection plugin.
	To be used for filenames of type tile_x{x}_y{y}.tif, in combination with stitching with known coordinates.
	See getFilenamePatternOfTiles() for alternative filenames
	Note: assumes that the counter is 1,2,3,... , i.e. not 01,02,03.. or 4,5,6,...
		and assumes that the filenames have exactly this structure: fffff_x1fffff_y1ffffff.tif etc. ("f" is placeholder (letter, or _, etc) but cannot be "_x","_y").
		Infers filename pattern from first file in input_dir.
	"""
	x_identifier="_x"
	y_identifier="_y"

	fullfiles=getIntegerSortedFileList(input_dir,extension)

	# infer name from first file in list
	basename=fullfiles[0].split(os.sep)[-1] # incl. extension

	# find prefix
	xloc=basename.find(x_identifier)
	prefix=basename[:xloc]+x_identifier # e.g. "img_patch_x"

	# find intermediate
	yloc=basename.find(y_identifier)
	middleix=basename[xloc+1+len(x_identifier):yloc]+y_identifier # e.g. "_y". assumes that the first filename has a 1-digit tile id!

	# find suffix
	suffix=basename[yloc+1+len(y_identifier):]# e.g. ".tif". assumes that the first filename has a 1-digit tile id!

	filename_pattern=prefix+"{x}"+middleix+"{y}"+suffix # grid collection syntax
	return filename_pattern



# ==================================================================================
# surface_extraction step: main functionality
# ==================================================================================
# Loads a z stack (typically full fused image or a roi) and finds the mincostzsurface. 
# Then reslices along this surface and saves resliced stack + heightmap.
# THe script reads all files which are in the output folder of stitching: -> 2b_stitched
# If a roi should be processed then save a roi of the fused image in that folder and adjust filename_must_contain below.
#

def extractFlattenedSurfaces(input_dir, save_dir, downsamplefactor_xy=0.25, downsamplefactor_z=1, max_dz=5, 
	numslices_above=5, numslices_below=5, sigma_cost=0, sigma_zmap=11, extension="tif", patch_increment=-1):
	"""
	Calculates the surface (zmap) of the stitched image and then extract flattened layers along that surface.
	Core plugin that is used: mincostzsurface
	Params:
	downsamplefactor_xy # mincostzsurface: downsampling for surface extraction. 
	downsamplefactor_z # mincostzsurface: downsampling for surface extraction. 
	max_dz: mincostzsurface: max change in z for neighboring pixels (creates smoothness constraint)
	numslices_above: number of slices of resliced stack: numslices_above+numslices+below=1 (the centerslice)
	numslices_below
	sigma_cost: (Deprecated).gaussian blur of the cost function (zstack) to extract a smoother surface (avoid muscle fiber pattern). (active if >0). 
	sigma_zmap: gaussian blur of the calculated height map to extract a smoother surface (avoid muscle fiber pattern). (active if >0) # NOTE: Applying 2 sigmas could be redundant ...
	patch_increment: from global variables
	"""
	global percentage_patch_overlap
	
	IJ.log("\nStarting extractFlattenedSurfaces ...")

	# prepare folders 
	assert os.path.exists(input_dir)
	if not os.path.exists(save_dir): 
		os.makedirs(save_dir)

	# save directories for substeps
	save_dir_imgpatches=os.path.join(save_dir,"image_patches") # input image patches
	save_dir_heightpatches=os.path.join(save_dir,"heightmap_patches") # retrieved heightmap per patch
	save_dir_surfacepatches=os.path.join(save_dir,"flattenedsurface_patches") # retrieved heightmap per patch

	if not os.path.exists(save_dir_imgpatches):
		os.makedirs(save_dir_imgpatches)
	if not os.path.exists(save_dir_heightpatches):
		os.makedirs(save_dir_heightpatches)
	if not os.path.exists(save_dir_surfacepatches):
		os.makedirs(save_dir_surfacepatches)

	# get stitched image
	filename_imgstack=os.path.join(input_dir,title_fused_image+"."+extension)
	imp=BF.openImagePlus(filename_imgstack)[0]

	# extract calibration
	calib=imp.getCalibration()

	# handle -1 case: use full image
	if patch_increment<0:
		patch_increment=10*max(imp.getWidth(),imp.getHeight());

	# resave in patches
	IJ.log("Resaving image into patches")
	nx,ny, patch_size=splitImageIntoPatchesAndSave(imp, save_dir_imgpatches, patch_increment,percentage_patch_overlap,verbose=True)

	# update the actual overlap (compensate for integer rounding)
	percentage_patch_overlap=100*(1-float(patch_increment)/patch_size)

	# compute heightmap per patch
	patchfiles=getIntegerSortedFileList(save_dir_imgpatches,"tif") # integer sorting!
	for idx in range(len(patchfiles)):
		# process image
		IJ.log("Computing height map of image: "+ patchfiles[idx])
		tstart=time.time()

		patch_imp=BF.openImagePlus(patchfiles[idx])[0]

		# compute height map (slow)
		computeHeightMapAndSave(patch_imp, save_dir_heightpatches, downsamplefactor_xy, downsamplefactor_z, max_dz, sigma_cost, sigma_zmap)
		
		IJ.log("Elapsed time:"+ str( int(time.time()-tstart) )+" sec.")

	# extract the surface per patch. Note: we could interleave this directly with the loop above
	heightmapfiles=getIntegerSortedFileList(save_dir_heightpatches,"tif") # integer sorted!
	assert len(heightmapfiles)==len(patchfiles)

	for idx in range(len(patchfiles)):
		patch_imp=BF.openImagePlus(patchfiles[idx])[0]
		patch_heightmap=BF.openImagePlus(heightmapfiles[idx])[0]
		resliceSurfaceAndSave(patch_imp,patch_heightmap,save_dir_surfacepatches,numslices_above,numslices_below)

	
	# stitch the extracted surface patches and save
	if len(patchfiles)==1:
		# nothing to stitch: copy the file
		surface_file=getIntegerSortedFileList(save_dir_surfacepatches,"tif")[0]
		surface_imp=BF.openImagePlus(surface_file)[0]

	else:
		filepattern=getFilenamePatternOfTilesXY(save_dir_surfacepatches,"tif")
		IJ.log("Stitching files with file pattern: "+str(filepattern))
		
		startidx=1

		# stitching: currently do compute overlap
		IJ.run(None, "Grid/Collection stitching", "type=[Filename defined position] order=[Defined by filename         ] grid_size_x="+str(nx)+" grid_size_y="+str(ny)+ \
			   " tile_overlap="+str(percentage_patch_overlap)+" first_file_index_x="+str(startidx)+" first_file_index_y="+str(startidx)+" directory="+save_dir_surfacepatches+ \
			   " file_names=["+filepattern+"] output_textfile_name=TileConfiguration.txt" +	" fusion_method=[Linear Blending] regression_threshold=0.3 max/avg_displacement_threshold=2.50"+ \
			   " absolute_displacement_threshold=3.50 " +\
			   " compute_overlap " +\
			   " computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
		# Note: in an earlier version the heightmaps wree stitched. This failed because the value of the height-patches was changed during the stitching (intensity adjustement?)
		
		# grab stitched image
		surface_imp=IJ.getImage()
		surface_imp.hide()

		# recover calibration
		surface_imp.setCalibration(calib)

	# save the fused surface image
	IJ.saveAsTiff(surface_imp,os.path.join(save_dir,"flattened_surface_of_full_zstack_stitched.tif"))

	IJ.log("Finished surface extraction.\n")



# ==================================================================================
# surface_extraction step: helper functions
# ==================================================================================
def correctImageSize2D(imp, targetWidth, targetHeight):
	"""padding/cropping along width & height to adjust the image dims to the target dims"""

	# ensure that new dims are less or equal to target dims
	roi=Roi(0,0,min(imp.getWidth(),targetWidth),min(imp.getHeight(),targetHeight))
	imp.setRoi(roi)
	imp=Duplicator().run(imp)

	# ensure that new dims are equal to target dims: enlarge
	newWidth=imp.getWidth()
	newHeight=imp.getHeight()

	paddingX=max(0,targetWidth-newWidth) # can be 0
	paddingY=max(0,targetHeight-newHeight)

	IJ.run(imp,"Extend Image Borders", "left=0 right="+str(paddingX)+" top=0 bottom="+str(paddingY)+" fill=Replicate");
	imp=IJ.getImage();
	imp.hide()

	return imp


def splitImageIntoPatchesAndSave(imp, save_dir, patch_increment=500, percentage_patch_overlap=15, verbose=False):
	"""splits an image imp into patches and saves them into save_dir. If the patch_increment does not fit into image size as an integer the patches of the last row & column can be larger.
	imp: 2D or 3D (x,y,z)
	patch_increment: increment of start coordinates of patches (along x and y, in pixel). This is the net size, the additional size is the overlap with neighboring patches.
	percentage_patch_overlap. between [0.100]. should be ca 10-20. Full patch size is: patch_increment/(1-percentagepatch_overlap/100)
	returns: nx,ny: the number of patches along x & y direction
			patch_size: the patch size in pixel including the overlap
	"""
	assert (percentage_patch_overlap>0 and percentage_patch_overlap<100)

	width,height,channels,slices,frames=imp.getDimensions();

	# the actual patch size incl the overlap
	patch_size=int(patch_increment/(1-float(percentage_patch_overlap)/100.0))

	# top-left corners of the patches
	xstart_list=[0]
	for x in range(patch_increment,width-patch_increment,patch_increment):
		xstart_list.append(x)

	ystart_list=[0]
	for y in range(patch_increment,height-patch_increment,patch_increment):
		ystart_list.append(y)

	nx=len(xstart_list)
	ny=len(ystart_list)

	IJ.log("Splitting the image into patches: nx="+str(len(xstart_list))+ ", ny="+str(len(ystart_list)))

	# width and height of the patches
	width_list=[]
	for idx in range(nx):
		if idx<nx-1:
			width_list.append(patch_size)
		else: # extend until the image border
			width_list.append(width-xstart_list[-1])

	height_list=[]
	for idx in range(ny):
		if idx<ny-1:
			height_list.append(patch_size)
		else: # extend until the image border
			height_list.append(height-ystart_list[-1])

	if verbose:
		IJ.log("xstart coordinates: "+str(xstart_list))
		IJ.log("ystart coordinates: "+str(ystart_list))
		IJ.log("patch widths:       "+str(width_list))
		IJ.log("patch heights:       "+str(height_list))


	# create the patches and save them
	for xidx in range(nx):
		for yidx in range(ny):
			imp.killRoi();
			roi=Roi(xstart_list[xidx],ystart_list[yidx],width_list[xidx],height_list[yidx])
			imp.setRoi(roi);
			cropImp=Duplicator().run(imp)
			imp.killRoi()

			# save with a proper tile name (for grid/collection stitching)
			fn="imgpatch_x"+str(1+xidx)+"_y"+str(1+yidx)+".tif"
			IJ.saveAsTiff(cropImp, os.path.join(save_dir,fn));

	return nx,ny, patch_size




def resliceSurfaceAndSave(imp, heightmap, save_dir, numslices_above, numslices_below):
	"""
	Reslices the image imp along the precomputed heightmap and saves the result. (part 2 of minCostZSurface)
	For parameters see global variables on top of script.
	"""
	title=imp.getTitle()

	# reslice
	IJ.log("Reslicing along surface")
	imp_surface = resliceStackAlongSurface(imp, heightmap, numslices_above, numslices_below)

	# if more than one channel: z & c are flipped
	channels=imp_surface.getNChannels()
	if channels>1:
		imp_surface=HyperStackConverter.toHyperStack(imp_surface,1,channels,1)

	IJ.save(imp_surface,os.path.join(save_dir,"flattened_surface_"+imp.getTitle().rsplit(".",1)[0]+".tif"))

	# close windows
	IJ.run(imp, "Close All", "")





def computeHeightMapAndSave(imp, save_dir, downsamplefactor_xy, downsamplefactor_z, max_dz, sigma_cost, sigma_zmap):
	"""
	Computes heightmap, part of surface reconstruction with mincostzsurface, then saves the computed map.
	Typically imp is either the full stitched image, or image patches (from splitImageIntoPatches)
	save_dir: should be a subfolder only for the heightmaps
	Height map and cost function are (can be) blurred for smoother reconstruction. For parameters see global variables on top of script.
	"""
	title=imp.getTitle()

	# create cost image
	imp_cost=makeCostImage(imp, sigma_cost)

	# calculate height map (imglib2 type)
	IJ.log("Creating Z surface (very slow)")
	imgzmap= createZmap(imp_cost,downsamplefactor_xy,downsamplefactor_z, max_dz)
	imp_zmap=ImageJFunctions.wrap(imgzmap,"height map")
	if sigma_zmap>0:
		IJ.run(imp_zmap, "Gaussian Blur...", "sigma="+str(sigma_zmap)+" stack")

	# reslice
	IJ.log("Computing height map")
	imp_surface = resliceStackAlongSurface(imp, imgzmap, numslices_above, numslices_below)

	# if more than one channel: z & c are flipped
	channels=imp_surface.getNChannels()
	if channels>1:
		imp_surface=HyperStackConverter.toHyperStack(imp_surface,1,channels,1)

	IJ.save(imp_zmap,os.path.join(save_dir,"heightmap_"+imp.getTitle().rsplit(".",1)[0]+".tif"))

	# close windows
	IJ.run(imp, "Close All", "")



def flattenAndSave(imp, save_dir, downsamplefactor_xy, downsamplefactor_z, max_dz, numslices_above, numslices_below, sigma_cost, sigma_zmap):
	""" DEPRECATED
	Surface reconstruction with mincostzsurface, then reslicing and saving.
	Height map and cost function are (can be) blurred for smoother reconstruction. For parameters see global variables on top of script.
	"""
	title=imp.getTitle()

	# create cost image 
	imp_cost=makeCostImage(imp, sigma_cost)

	# calculate height map (imglib2 type)
	IJ.log("Creating Z surface (very slow)")
	imgzmap= createZmap(imp_cost,downsamplefactor_xy,downsamplefactor_z, max_dz)
	imp_zmap=ImageJFunctions.wrap(imgzmap,"height map")
	if sigma_zmap>0:
		IJ.run(imp_zmap, "Gaussian Blur...", "sigma="+str(sigma_zmap)+" stack")

	# reslice
	IJ.log("Reslicing along surface")
	imp_surface = resliceStackAlongSurface(imp, imgzmap, numslices_above, numslices_below)

	# if more than one channel: z & c are flipped
	channels=imp_surface.getNChannels()
	if channels>1:
		imp_surface=HyperStackConverter.toHyperStack(imp_surface,1,channels,1)
		
	# save results
	save_dir_heightmaps=os.path.join(save_dir,"heightmaps")
	if not os.path.exists(save_dir_heightmaps):
		os.makedirs(save_dir_heightmaps)

	IJ.save(imp_surface,os.path.join(save_dir,"flattened_surface_"+imp.getTitle().rsplit(".",1)[0]+".tif"))
	IJ.save(imp_zmap,os.path.join(save_dir_heightmaps,"heightmap_"+imp.getTitle().rsplit(".",1)[0]+".tif"))

	# close windows
	IJ.run(imp, "Close All", "")
		


def makeCostImage(imp, sigma=0):
	"""Creates a cost image by duplicating + inverting + optional smoothing.
	data type: imageplus
	"""
	cost_imp = imp.duplicate()

	IJ.run(cost_imp, "Invert", "stack")
	if sigma>0:
		IJ.run(cost_imp, "Gaussian Blur...", "sigma="+str(sigma)+" stack")

	return cost_imp
	

def createZmap(cost_imp,downsample_factor_xy=1, downsample_factor_z=1, max_dz=1):
	"""Calls MinCostZSurface (1 surface).
	cost_imp: imageplus
	zMap: imglib2 image
	"""	
	# convert to imglib2 img 
	cost_img = ImageJFunctions.wrap(cost_imp)

	zMap = op.run( "MinCostZSurface",
				cost_img			,
				downsample_factor_xy,
				downsample_factor_z	,
				max_dz 	            )

	return zMap

	
def resliceStackAlongSurface(imp, imgzmap, slices_above, slices_below):
	"""
	imp: input image. imageplus type
	imgzmap: height map. imglib2 type
	returns: imp_surface:
	"""
	Img= ImageJFunctions.wrap(imp)
	ImgSurface = op.run( "zMapReslice", Img, imgzmap , slices_above ,slices_below)

	# convert to imageplus
	imp_surface=ImageJFunctions.wrap(ImgSurface,"flattened surface")
	
	return imp_surface
	

# ==================================================================================
# derotation step: main functionality
# ==================================================================================
# derotates the fused image and the zmax projection of the fused image

def derotateImages(input_dir_stitch, input_dir_surface, save_dir, derotation_mode, zmax_identifier="zmax",flatsurface_identifier="flatten",extension="tif"):
	"""Derotates images in input_dir_stitch and input_dir_surface such that flatworm gets vertically aligned (head at the top).
	ZMax projected image (from input_dir_stitch) is used for Rotation determination.
	These images are then derotated: stitched stack (if existent), stitched zmax, extracted surface (also stitched and used for further processing.)
	Params:
	derotation_mode: 0: no derotation, 1: automatic axis calculation + derotation (very basic!!), 2: interactive by user
	flatsurface_identifier: string to distinguish between flattened surface and heightmap
	zmax_identifier: (sub)string which is unique to the zmax image in the input_folder
	"""
	IJ.log("\nStarting derotateImages ...")
	
	# prepare folders 
	assert os.path.exists(input_dir_stitch)
	assert os.path.exists(input_dir_surface)
	if not os.path.exists(save_dir):
		os.makedirs(save_dir)

	# get worm rotation (deviation from horizontal line)
	if derotation_mode==0:
		IJ.log("Rotating by 90deg. Skipping fine adjustement of derotation (derotation_mode=0).")
		angle=None
	elif derotation_mode==1:			
		IJ.log("Automatic detection of rotation. Simple method(!).")
		angle=extractRotationAutomatically(input_dir_stitch, zmax_identifier)
		# avoid 'flip' rotations of ca 180deg (extracted angle from fitEllipse has +-180deg ambiguity. keeping the angle at [-90deg,90deg] is best possible rotation guess)
		if angle is not None and angle>90:
			angle=angle-180
			IJ.log("corrected derotation angle: "+str(angle))
	else: 
		IJ.log("Interactive derotation")
		angle=extractRotationInteractively(input_dir_stitch, zmax_identifier) # angle between -180 and 180
	
	# apply derotation to all images in folder
	# find fused image and zmax_image
	imagefiles=getFileList(input_dir_stitch,extension=extension)
	# find the flattened_surface image
	imagefilestmp=getFileList(input_dir_surface,extension=extension)
	for filename in imagefilestmp:
		basename=filename.split(os.sep)[-1]
		if flatsurface_identifier is None or flatsurface_identifier in basename:
			imagefiles.append(filename)


	for idx in range(len(imagefiles)):
		# get the basename
		filename=imagefiles[idx]
		filename=filename.split(os.sep)[-1]
		filename=filename.split(".")[0] # rm extension

		# open image
		IJ.log("\nProcessing image "+ imagefiles[idx])
		
		imp=BF.openImagePlus(imagefiles[idx])[0] 
		imp.killRoi() # image can only be enlarged if the line roi is not visible (otherwise attempts to rotate roi only)

		# derotate
		imp.show() # bug in Rotator: stack rotation doesn't work (only half of the stack rotated) if imp not visible

		# ... by 90deg for vertical alignment
		IJ.run(imp, "Rotate 90 Degrees Left", "");

		# .. . by deviation from horizontal axis
		if angle is not None:
			IJ.run(imp, "Rotate... ", "angle="+str(angle)+" grid=1 interpolation=Bilinear enlarge stack") # clockwise
		
		# save
		IJ.save(imp,os.path.join(save_dir,filename+"_derotated.tif"))

		imp.changes=False
		imp.close() 
			

	IJ.log("Finished derotation.\n")


# ==================================================================================
# derotation step: helper functions
# ==================================================================================
def extractRotationAutomatically(input_dir,zmax_identifier):
	""" Rotation based on zmax image. Loads zmax image, thresholds + morph ops,
	determines angle of major/minor axis ellipse -> rotation angle.
	Simple algorithm.
	"""
	
	# find & open image with zmax projection
	try:
		fn_zmax=glob.glob(input_dir+os.sep+"*" +zmax_identifier+"*")[0]
	except:
		IJ.log("No ZMax image found. Skipping derotation.")
		return None
		
	imp=BF.openImagePlus(fn_zmax)[0] 
	
	# get a simple mask of the (interior) flatworm region
	IJ.run(imp, "Gaussian Blur...", "sigma=5")
	#IJ.setAutoThreshold(imp, "Triangle dark no-reset")
	#IJ.run(imp, "Convert to Mask", "");
	IJ.run(imp, "Auto Threshold", "method=Triangle ignore_black white")  # via different threshold window!
	IJ.run(imp, "Options...", "iterations=10 count=1 black do=Dilate") # connect fiber masks
	IJ.run(imp, "Options...", "iterations=1 count=1 black do=Nothing") # reset iterations
	IJ.run(imp, "Fill Holes", "");
	IJ.run(imp, "Options...", "iterations=20 count=1 black do=Erode")
	IJ.run(imp, "Options...", "iterations=1 count=1 black do=Nothing") # reset iterations
	IJ.run(imp, "Create Selection", "")

	stats=imp.getStatistics(Measurements.ELLIPSE)
	IJ.log("Detected angle (automatic mode): "+str(stats.angle))
	return stats.angle



def extractRotationInteractively(input_dir,zmax_identifier):
	""" Rotation based on zmax image. Loads zmax image, then asks user to draw a line along the long axis.
	Rotation is the angle of the line.
	"""
	
	# find & open image with zmax projection
	try:
		fn_zmax=glob.glob(input_dir+os.sep+"*" +zmax_identifier+"*")[0]
	except:
		IJ.log("No ZMax image found. Skipping derotation.")
		return None
		
	imp=BF.openImagePlus(fn_zmax)[0] 
	imp.show()

	WaitForUserDialog("User input required","Draw a straight line along the centerline of the worm (from tail to head).\nThen click ok.").show()
	roi=imp.getRoi()
	imp.close()

	if roi is not None and roi.getType()==Roi.LINE:
		angle=roi.getAngle() 
		IJ.log("Detected angle (interactive mode): "+str(angle))
		return angle
	else:
		IJ.log("A Straight Line selection is required for derotation. Skipping derotation.")
		return None



# ==================================================================================
# orientation_analysis step: main functionality
# ==================================================================================
# Does Orientation analysis on a surface-flattened z-stack (core plugin: OrientationJ). 
# Loads a flattened z stack and does orientation analysis per slice. saves color coded image and histograms.
# Includes functionality for simple masking to obtain histogram on roi region.
# The script reads all files which are in the output folder of surface_extraction: -> 4_derotation
# If a roi should be processed then save a roi of the fused image in that folder and adjust filename_must_contain below.
#
def analyzeFiberOrientation(input_dir, save_dir, filename_must_contain, filename_must_not_contain, mask_mode, vector_grid_size=50, extension="tif"):
	"""
	Does Orientation analysis on a surface-flattened z-stack (core plugin: OrientationJ). 
	Loads a flattened z stack and does orientation analysis per slice. saves color coded image and histograms.
	Includes functionality for simple masking to obtain histogram on roi region.
	If a roi should be processed then save a roi of the fused image in that folder and adjust filename_must_contain below.
	Core plugin that is used: OrientationJ
	Update 2019-03-12:
	* The orientation is now also displayed as vector field. Note that the vector field plugin from OrientationJ had
		bugs (no macro access) so it was in  parts re-coded here.
		Check regularly the OrientationJ code (or update site) to check if the bug gets fixed.
		Length of the vectors in vectorField: propto coherency
	*  Histogram now weighed with coherency
	Params:
	filename_must_contain: should be "flattened_surface". if None: no selection
	filename_must_not_contain: should be "zmax". if None: no selection
	mask_mode: 0: no masking (use full image), 1: autodetermined (simple thresholding mask), 2: drawn by user
				determines region in which orientation histogram is calculated
	vector_grid_size:  px. grid size (length & width) of patches for orientation averaging for vectorfield plot
	

	Update 2019-03: the results in raw format (orientation image, roi, resutls tables) are now stored in a subfolder /data/
		to allow for further processing in python etc.
	"""
	IJ.log("\nStarting analyzeFiberOrientation ...")
	
	# prepare folders 
	assert os.path.exists(input_dir)
	if not os.path.exists(save_dir): 
		os.makedirs(save_dir)

	# folder for storing results data in raw format
	data_savedir=os.path.join(save_dir,"data")
	if not os.path.exists(data_savedir):
		os.makedirs(data_savedir)


	# get list with all imagefiles in input folder (typically 1 image)
	imagefiles=getFileList(input_dir,extension=extension, verbose=False)
	
	# loop over all files
	for idx in range(len(imagefiles)):
		# get the basename
		filename=imagefiles[idx]
		filename=filename.split(os.sep)[-1]

		# process image
		if (filename_must_contain is None) or (filename_must_contain in filename):
			if (filename_must_not_contain is None) or not (filename_must_not_contain in filename):
				IJ.log("\nProcessing image "+ imagefiles[idx])
				tstart=time.time()

				bfoptions=ImporterOptions()
				bfoptions.setId(imagefiles[idx])
				imp=BF.openImagePlus(bfoptions)[0]
				assert imp.getNChannels()==1,"Image is multi-channel or other (e.g. default) channel order!"

				window_size=2.0 # window size of structure tensor (local orientation)

				# main step 1: draw roi, colored orientation image, histograms
				regionroi=analyzeOrientationAndSave(imp, save_dir, data_savedir, mask_mode, window_size) # returns the regionroi to be reused in main step 2 (vectorfield)

				# main step 2: vector field
				visualizeVectorFieldAndSave(imp, save_dir, data_savedir, regionroi, vector_grid_size, window_size)

				# close windows
				# IJ.run(imp, "Close All", "") # bug work around: allow for manual saving


				IJ.log("Elapsed time:"+str( int(time.time()-tstart))+" sec.")

			else:
				pass
				#IJ.log("Skipping image (does not contain required substring):"+ imagefiles[idx])

	IJ.log("Finished orientation analysis.\n")



# ==================================================================================
# orientation_analysis step: helper functions
# ==================================================================================
def analyzeOrientationAndSave(imp, save_dir, data_savedir, mask_mode, window_size):
	"""
	Runs "OrientationJ Analysis" analysis plugin and saves the orientation-color images, plus calculates histogram in
	a roi region (interactive/manual depending on mask_mode). Histogram is weighted by coherency.
	Saves results (RGB, flattened, contrast enhanced)
	:param imp: stack with flattened surfaces
	:param save_dir: for viz images, overlays, ...
	:param data_savedir: for tables, masks, orientation & coherency image. everthing for further processing in different software
	:param mask_mode: 0,1,2
	:param window_size: size of local window ('tensor') for orientation determination
	:return: roi: the roi of planarian, or None
	"""
	title=imp.getTitle()
	calib=imp.getCalibration();

	# get roi for histogram processing or set roi to none (do before orientaiton analysis to keep the ointeractive steps close together)
	roi, imp_proj=getRoiForHistogram(imp, mask_mode)
	
	# == Orientation analysis ==
	# OrientationJ provides automatic Macro scripting but not python scripting. 
	# Therefore, wrap the macro to execute it here. Requires processed image to be visible.
	
	# result window names (defaults of OrientationJ)
	titlecolor="OJ-Color-survey-1" 
	titleorientation="OJ-Orientation-1" 
	titlecoherency="OJ-Coherency-1" 
	
	# bring correct window for processing to the front (macro requires image to be visble)
	imp.show()
	win=WindowManager.getWindow(title)
	WindowManager.setCurrentWindow(win)
	
	# == run OrientationJ_Analysis plugin ==
	IJ.log("starting OrientationJ Analysis (slow)")
	mr=Macro_Runner()
	# notes:  gradient type 0: cubic spline interpol. orientation window must be in degree (not rad)
	mr.runMacro("run(\"OrientationJ Analysis\", \"tensor="+str(window_size)+" gradient=0 color-survey=on hsb=on hue=Orientation sat=Coherency bri=Original-Image orientation=on coherency=on radian=off \");","")
	IJ.log("Compututation of orientation and coherency map finished")
	
	# == grab color image, rescale intensity and save ==
	IJ.log("Creating color image")
	imp_color=WindowManager.getImage(titlecolor)
	imp_color.setCalibration(calib)
	
	saturated_percentage=0.35
	if imp_color.getNSlices()==1: # single image
		IJ.run(imp_color, "Enhance Contrast...", "saturated="+str(saturated_percentage))
	else: # use stack hist + processing
		IJ.run(imp_color, "Enhance Contrast...", "saturated="+str(saturated_percentage)+" process_all use");
	IJ.save(imp_color,os.path.join(save_dir,"coloredOrientation_RGB_contrastEnhanced.tif"))
	imp_color.close()
	
	# == save more images ==
	# save a mask image of the roi (for downstream processing in other software)
	mask=createMask(imp,roi)
	IJ.save(mask,os.path.join(data_savedir,"roiMask.tif"))

	# save an overlay of the mask outline (for visualization)
	IJ.run(imp_proj, "Enhance Contrast", "saturated=0.35")
	IJ.run(imp_proj, "8-bit", "")
	if roi is not None:
		ov=Overlay()
		ov.add(roi)
		imp_proj.setOverlay(ov)
	imp_proj=imp_proj.flatten()
	IJ.save(imp_proj,os.path.join(save_dir,"roiOutline_RGB_contrastEnhanced_zmax_of_flattenedsurface.tif"))

	# save the orientation & coherency images in 32bit for potential downstream processing
	imp_coh=WindowManager.getImage(titlecoherency)
	IJ.save(imp_coh,os.path.join(data_savedir,"Coherency.tif"))
	imp_ori=WindowManager.getImage(titleorientation)
	IJ.save(imp_ori,os.path.join(data_savedir,"Orientation_degrees.tif"))

	# == get orientation histogram in a roi region and save ==
	# grab images
	IJ.log("creating weighted histogram of orientations (very slow)")
	imp_orientation=WindowManager.getImage(titleorientation)
	imp_coherency=WindowManager.getImage(titlecoherency)
	IJ.run(imp_orientation,"32-bit","")

	# tables for collecting the histogram distributions
	rt_weighted=ResultsTable() # recommended
	rt_unweighted=ResultsTable()
	rt_weighted.setPrecision(10)
	rt_unweighted.setPrecision(10)

	# collect histogram (once coherency weighted, once without weighting)
	for sliceidx in range(1,1+imp_orientation.getNSlices()):
		imp_orientation.killRoi()
		imp_coherency.killRoi()
		impOslice=Duplicator().run(imp_orientation,sliceidx,sliceidx)
		impCslice=Duplicator().run(imp_coherency,sliceidx,sliceidx)

		# coherency-weighted orientation (recommended for further analysis)
		hist_weighted=calculateHistogram(impOslice,mask,impCslice,normalize=True)
		fig=createHistPlot(hist_weighted, True, sliceidx)
		#IJ.save(fig,os.path.join(save_dir,"Histogram_CoherencyWeighted_slice"+str(sliceidx)+".png"))
		IJ.save(fig,os.path.join(save_dir,"Histogram_slice"+str(sliceidx)+"CoherencyWeighted.png"))
		 
		# unweighted orientation histogram.  
		hist_unweighted=calculateHistogram(impOslice,mask,None, normalize=True)
		fig=createHistPlot(hist_unweighted, False, sliceidx)
		#IJ.save(fig,os.path.join(save_dir,"Histogram_UnWeighted_slice"+str(sliceidx)+".png"))
		IJ.save(fig,os.path.join(save_dir,"Histogram_slice"+str(sliceidx)+"_Unweighted.png"))
		
		# collect hist's in results table
		rt_weighted=addSliceHistToTable(rt_weighted,hist_weighted, sliceidx)
		rt_unweighted=addSliceHistToTable(rt_unweighted,hist_unweighted, sliceidx)
		
	# save histogram tables
	rt_weighted.saveAs(os.path.join(data_savedir,"orientationDistribution_CoherencyWeighted.csv"))
	rt_unweighted.saveAs(os.path.join(data_savedir,"orientationDistribution_Unweighted.csv"))

	imp_orientation.changes=False # suppress pop-up
	imp_orientation.close()
	imp_coherency.close()

	return roi # for vector field visualization


def getRoiForHistogram(imp, mask_mode):
	"""	 get roi for histogram processing or set roi to none"""
	# get a 2D image version
	if imp.getNSlices()>1:
		imp_proj=ZProjector.run(imp,"max") # a projection of a zmapresliced flattened stack! for segmentation only. no biological meaning
	else: # for development
		imp_proj=imp.duplicate()	

	if mask_mode==0:
		roi=None
	elif mask_mode==1:
		IJ.log("Automatic creation of roi mask")
		# get a simple mask of the (interior) flatworm region
		maskr=imp_proj.duplicate()
		IJ.run(maskr, "Gaussian Blur...", "sigma=5")
		IJ.run(maskr, "Auto Threshold", "method=Triangle ignore_black white")  # via different threshold window!
		IJ.run(maskr, "Options...", "iterations=10 count=1 black do=Dilate") # connect fiber masks
		IJ.run(maskr, "Options...", "iterations=1 count=1 black do=Nothing") # reset iterations
		IJ.run(maskr, "Fill Holes", "");
		IJ.run(maskr, "Options...", "iterations=10 count=1 black do=Erode")
		IJ.run(maskr, "Options...", "iterations=1 count=1 black do=Nothing") # reset iterations
		IJ.run(maskr, "Create Selection", "")
		roi=maskr.getRoi()
		maskr.killRoi()
	else:
		IJ.log("Interactive creation of roi mask")
		imp_proj.show()
		WaitForUserDialog("User input required","Draw a ROI on the MAX-projected image (region for orientation histogram calculation). Then click ok.").show()
		roi=imp_proj.getRoi()
		imp_proj.killRoi()
		if roi is None:
			IJ.log("No roi was specified. Will process full image.")

	return roi, imp_proj


def visualizeVectorFieldAndSave(imp, save_dir, data_savedir, regionroi, gridSize=50, windowSize=2.0):
	"""
	Runs "OrientationJ Vector Field" to create a nematic vector field of the fiber orientations.
	The core orientation algorithm has to be rerun (after analyzeOrientationAndSave) to get access to the other output
	variables (not very efficient ...)
	Vector length is proportional to orientation coherency.
	Results are saved to disk (RGB, flattened, contrast enhanced).

	IMPORTANT: The VectorField plugin has bugs: It could not be run as macro. Therefore, we access OrientationJ at a lower level here.
			It also had errors in patch averaging (for the nematic vector calculation), this was filed as github issue and is fixed now.

	:param imp: stack with flattened surfaces
	:param save_dir: for viz images, overlays, ...
	:param: data_savedir: for tables, masks, ...
	:param regionroi: from analyzeOrienationAndSave. a Roi or None
	:param gridSize:  px. grid size (length & width) of patches for orientation averaging
	:param windowSize:  also called 'tensor' argument in the orientation GUI. size of local region for orientation calculation. must match analysis plugin
	"""

	imp.show()
	title=imp.getTitle()

	vectorScale=90.0 # %. max vector length as percentage of gridSize
	strokeColor=Color.cyan
	strokeWidth=3

	# run the vector field plugin HACK VERSION. it acts on the current image
	IJ.log("Computing Vector field (very slow).")
	# #run_VectorFieldPlugin_MediumLevelAccess(windowSize, gridSize, vectorScale) # used until v1.0. caused mysterious single slice vectorfield?
	run_VectorFieldPlugin_LowLevelAccess(windowSize, gridSize, vectorScale)

	# grab the results table with the vector field data
	tableTitle="OJ-Table-Vector-Field-"
	rt = WindowManager.getWindow(tableTitle).getTextPanel().getResultsTable()
	time.sleep(3)

	# save the rt for re-use in python (beware of the flipped -dx!)
	rt.saveAs(os.path.join(data_savedir,"vectorFieldData.csv"))

	# draw all vectors within the roi (as overlay)
	drawRestrictedVectorFieldOverlay(imp, regionroi,gridSize, rt=rt, strokeColor=strokeColor,strokeWidth=strokeWidth)

	# close results table
	# WindowManager.getWindow(tableTitle).close() # bugworkaround: allow for manual saving

	# add the user-drawn roi to the overlay
	if regionroi is not None:
		ov=imp.getOverlay()
		if ov is None:
			ov=Overlay()
		regionroi.setStrokeColor(Color.yellow)
		regionroi.setStrokeWidth(1)
		ov.add(regionroi);
		imp.setOverlay(ov)

	imp.killRoi()

	# create a flattened (imprinted) version for saving
	impvector_print=imp.duplicate()
	saturated_percentage=0.35
	if impvector_print.getNSlices==1: # single image
		IJ.run(impvector_print, "Enhance Contrast...", "saturated="+str(saturated_percentage))
	else: # use stack hist + processing
		IJ.run(impvector_print, "Enhance Contrast...", "saturated="+str(saturated_percentage)+" process_all use");
	IJ.run(impvector_print, "8-bit", "")

	if impvector_print.getNSlices()==1:
		impvector_print=impvector_print.flatten()
	else:
		impvector_print.flattenStack()

	# display iamges to avoid bug during saving (in windows7: otherwise only 1 overlay slice is saved)
	impvector_print.show()
	imp.show()

	# save images
	IJ.save(impvector_print,os.path.join(save_dir,"vectorFieldOrientation_RGB_contrastEnhanced.tif"))
	IJ.save(imp,os.path.join(save_dir,"vectorFieldOrientation_rawWithOverlay.tif"))





def calculateHistogram(impOrientation,mask=None,impCoherency=None, normalize=False):
	"""Calculates histogram of orientations, based on impOrientation.
	Histogram is from -90 to 90 with binSize=1 (-> 180bins)
	impOrientation: The orientation image (values must be in degree, not rad), output of OrientationJ analysis
	mask: only pixels with value>0 (=255) contribute. For None: all px contribute. For no roi: use white mask or None mask.
	If impCoherency is not None, then the histogram is a weighted histogram. Each pixel is weighted with its coherency.
			Values of coherency range from 0 to 1. impCoherency: output of OrientationJ analysis.
			It is recommended to always weigh with coherency
	normalize: if True then the histogram is normalized by total nr of pixels (of the roi).
		For unweighted histogram: bins sum then to =1
		For weighted histogram: bins sum to avg coherency within the roi (between 0 and 1)
	impOrientation, impCoherency: 2D images, no stack
	"""
	
	# prepare image	
	impOrientation.killRoi()
	IJ.run(impOrientation,"32-bit","")
	
	width=impOrientation.getWidth()
	height=impOrientation.getHeight()
	
	ipOrient=impOrientation.getProcessor()
	
	# histogram bins (binning interval hardcoded to be =1! (faster))
	histStart=-90
	histEnd=90
	histBins=range(histStart,histEnd,1) # start=lower values of bins (i.e. bin[0]: -90 to -89, bin[180]: 89 to 90
	
	# init the histogram count
	histCount=[0.0] * len(histBins)

	if impCoherency is not None:
		ipCoh=impCoherency.getProcessor()
		do_weighting=True
	else:
		do_weighting=False

	# needed for optional normalization
	pxcount=0

	if mask is not None:
		maskp=mask.getProcessor()
	
	# loop over all image pixels
	for x in range(width):
		for y in range(height):

			if mask is None or maskp.getPixelValue(x,y)>0:
				# orientation of this pixel (in deg!)
				angle=ipOrient.getPixelValue(x,y)
			
				# (weighted) contribution of this pixel
				if do_weighting:
					weight=ipCoh.getPixelValue(x,y)
				else:
					weight=1
		
				# corresponding index of bin in histogram
				histIdx=int(math.floor(angle-histStart))
				histIdx=min(len(histBins)-1,(max(histIdx,0)))
	
				# collect histogram
				histCount[histIdx]+=weight

				# for normalization
				pxcount+=1

	if normalize:
		histCount=[elem/float(pxcount) for elem in histCount]

	return histCount



def addSliceHistToTable(rt,histogram, sliceidx, angles=range(-90,90,1)):
			# angles must match the histogram creation angles!!!
			rt.incrementCounter()
			rt.addLabel("slice "+str(sliceidx))
			for idx in range(len(histogram)):
				rt.addValue(str(angles[idx]),histogram[idx])
			return rt


# @deprecated
def run_VectorFieldPlugin_MediumLevelAccess(windowSize=2.0, gridSize=50, vectorScale=90.0):
	"""	Recorded Macros of VectorField OrientationJ don't work (2019-03-12). 
	The issue is that https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/OrientationJ_Vector_Field.java
	does not have a "OrientationResults.show(process.getGroupImage(), params, 1);" at the end of the code.

	This function here mimicks the high-level plugin.

	In future: This recorded macro should be run (adjust the variable values)
	run("OrientationJ Vector Field", "tensor=2.0 gradient=0 radian=on vectorgrid=50 vectorscale=80.0 vectortype=2 vectoroverlay=off vectortable=on ");

	"""
	# see: https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/OrientationJ_Vector_Field.java
	from gui_orientation import WalkBarOrientationJ;
	from orientation import GroupImage;
	from orientation import OrientationParameters;
	from orientation import OrientationProcess
	from orientation import OrientationService;
	from orientation import OrientationResults # NEW

	params = OrientationParameters(OrientationService.VECTORFIELD);

	# the string parameters from the recording.
	# vectortype=2: scale length with coherency
	# windowSize: ='tensor' same as in OrientationJ Analysis
	macro_string_params="tensor="+str(windowSize)+" gradient=0 radian=off vectorgrid="+str(gridSize)+" vectorscale="+str(vectorScale)+" vectortype=2 vectoroverlay=off vectortable=on "

	params.getMacroParameters(macro_string_params);
	source = GroupImage.getCurrentImage();
	walk = WalkBarOrientationJ();
	process =  OrientationProcess(walk, source, params);
	process.run();

	# == show the results (missing in the original source code)
	OrientationResults.show(process.getGroupImage(), params, 1); # NEW



def run_VectorFieldPlugin_LowLevelAccess(windowSize=2.0, gridSize=50, vectorScale=90.0):
	""" This function replaces Workaround_to_run_vectorfield_orientationj_plugin() and accesses the same functionality at a lower level (see that function for more details on why it exists).
	This became necessary because for unknown reasons (large data size?) processing of stacks could fail: the vectorfield and table with vector stats was only displayed for the first slice.
	Here: * high level plugin is mimicked here to get access the plugin results 
		  * additionally: vectorField is calculated manually and explicitely for every single slice.
	returns: groupimage, orientationparameters
	"""
	# see: https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/OrientationJ_Vector_Field.java
	from gui_orientation import WalkBarOrientationJ;
	from orientation import GroupImage;
	from orientation import OrientationParameters;
	from orientation import OrientationProcess
	from orientation import OrientationService;
	from orientation import OrientationResults # NEW

	params = OrientationParameters(OrientationService.VECTORFIELD);

	# the string parameters from the recording.
	# vectortype=2: scale length with coherency
	# windowSize: ='tensor' same as in OrientationJ Analysis
	macro_string_params="tensor="+str(windowSize)+" gradient=0 radian=off vectorgrid="+str(gridSize)+" vectorscale="+str(vectorScale)+" vectortype=2 vectoroverlay=off vectortable=on "

	params.getMacroParameters(macro_string_params);
	source = GroupImage.getCurrentImage();
	walk = WalkBarOrientationJ();
	process =  OrientationProcess(walk, source, params);
	process.run();

	# == from here on: not in original OrientationJ code
	# manually compute the vectorfield
	table=reimplemented_computeVectorFieldForDisplay(process.getGroupImage(),params)

	# display the table. we could return it directly of course, but with displaying it is easier to switch between old and new code
	table.show("OJ-Table-Vector-Field-"); # now return directly



def reimplemented_computeVectorFieldForDisplay(gim, params):
	"""  reimplements the displayVectorField from here: https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/orientation/OrientationResults.java.
	reason: The original plugin sometimes only processes the first slice. it is not known whether the orientation analysis fails already or only the display vector field.
			here we start with re-implementing the display vectorfield part.
	NOTE: pretty slow since ported to python from java!
	gim: GroupImage
	params: OrientationParameters()

	returns: resultstable (i.e. table is not displayed but returned as variable)
	"""

	from orientation import Clusters
	from orientation import Cluster
	from java.lang import Math


	size = params.vectorGrid;
	type = params.vectorType;
	scale = params.vectorScale;

	gim.energy.getSizeZ();
	clusters = [Clusters() for i in range(gim.nt)]# syntax change from java->python

	xstart = ((gim.nx - (gim.nx / size) * size) / 2)
	ystart = int ((gim.ny - (gim.ny / size) * size) / 2)
	emax = gim.energy.getMaximum();
	if emax <= 0:
		return

	size2 = size * size;
	for t in range(gim.nt):
		clusters[t] = Clusters();
		for y in range(ystart, gim.ny, size):
			for x in range(xstart, gim.nx, size):
				dx=0.0
				dy=0.0
				coherencies=0.0
				energies=0.0

				for k in range(size):
					for l in range(size):
						angle = gim.orientation.getPixel(x+k, y+l, t);
						coh = gim.coherency.getPixel(x+k, y+l, t);
						dx += Math.cos(angle);
						dy += Math.sin(angle);
						coherencies += coh;
						energies += gim.energy.getPixel(x+k, y+l,t);

				dx /= size2;
				dy /= size2;
				coherencies /= size2;
				energies /= size2;
				if (energies > 0):
					if (coherencies > 0):
						clusters[t].add(Cluster(x, y, size, size, dx, dy, coherencies, (energies / emax)))

	# Table
	#if (params.showVectorTable): # now always calculate, and return directly
	table = ResultsTable();
	for t in range(gim.nt):
		for c in clusters[t]:
			a = Math.toDegrees(Math.atan2(c.dy, c.dx));
			if (a < -90):
				a += 180;
			if (a > 90):
				a -= 180;
			table.incrementCounter();
			table.addValue("X", c.x + size / 2);
			table.addValue("Y", c.y + size / 2);
			table.addValue("Slice", t);
			table.addValue("DX", -c.dx);
			table.addValue("DY", c.dy);
			table.addValue("Orientation", a);
			table.addValue("Coherency", c.coherency);
			table.addValue("Energy", c.energy);

	#table.show("OJ-Table-Vector-Field-"); # now return directly

	# Overlay: moved to separate function (with Roi filtering)

	return table






def drawRestrictedVectorFieldOverlay(imp, regionroi, gridSize, rt, vectorScale=90.0,strokeColor=Color.cyan, strokeWidth=2):
	"""Takes the results table of VectorField OrientationJ plugin and draws the vector lines.
	If a roi is given (e.g. planarian outline), drawing is restricted to the roi region.
	The results table must be open!
	Vector length is proportional to coherency

	params:
	imp: typically stack of flattened surfaces. can also be single slice. same image for which orientation plugin is run
	regionroi: area region or None
	gridSize: must be the same as for running the plugin. essential to compensate for shifts in center coord
	vectorScale: % of max vector length of grid size
	rt: ResultsTable with vectorfield

	"""
	# center points of vectors
	xCtrs=rt.getColumn(rt.getColumnIndex("X"))
	yCtrs=rt.getColumn(rt.getColumnIndex("Y"))

	# slice (z)
	zSlices=rt.getColumn(rt.getColumnIndex("Slice"))

	# coherency
	coherencies=rt.getColumn(rt.getColumnIndex("Coherency"))

	# angle (dx and dy)
	dxs=rt.getColumn(rt.getColumnIndex("DX"))
	dys=rt.getColumn(rt.getColumnIndex("DY"))

	# convert elements to int whenever known, and all other arrays to lists (to avoid surprises)
	xCtrs=[int(elem) for elem in xCtrs] # float->int, java array->list
	yCtrs=[int(elem) for elem in yCtrs]
	zSlices=[int(elem) for elem in zSlices]
	coherencies=list(coherencies)
	dxs=[-elem for elem in dxs] # dx was inverted when writing to the table! line 355 in https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/orientation/OrientationResults.java
	dys=list(dys)

	# create overlay. overlay creation based on original OrientationJ implementation:
	# https://github.com/Biomedical-Imaging-Group/OrientationJ/blob/master/src/main/java/orientation/OrientationResults.java
	# key modifications: the lines are displayed in the square center (instead of top-left corner), only vectors within a user-drawn roi are displayed, differently colored overlay
	overlay = Overlay()

	maxHalfLength=vectorScale/100.0 * gridSize * 0.5 # in px. of displayed lines

	for idx in range(len(zSlices)):

		# check if vector line is within user-roi
		if regionroi is None:
			pointInRoi=True
		else:
			pointInRoi=regionroi.contains(xCtrs[idx],yCtrs[idx])

		# only draw vector lines within user-roi
		if pointInRoi:
			l=maxHalfLength*coherencies[idx] # length proportional to coherency

			# end points of vector line
			x1=int(round(xCtrs[idx]+l*dxs[idx]))
			y1=int(round(yCtrs[idx]-l*dys[idx]))
			x2=int(round(xCtrs[idx]-l*dxs[idx]))
			y2=int(round(yCtrs[idx]+l*dys[idx]))

			# add line to overlay
			roi=Line(x1, y1, x2, y2)
			roi.setPosition(zSlices[idx] + 1);
			roi.setStrokeColor(strokeColor)
			roi.setStrokeWidth(strokeWidth)
			overlay.add(roi);
	
	imp.setOverlay(overlay)





def createHistPlot(hist,weighted,sliceidx):
	""" hist: list with (weighted) histogram count. bins must correspond to angle range(-90,90,1)
	weighted: True/False: determines title string
	sliceidx: integer: added to title string
	"""
	angles=[x for x in range(-90,90,1)]

	if weighted==True:
		figtitle = "OrientationHistogram_Masked_CoherencyWeighted_slice_"+str(sliceidx)
		color=Color.blue
	else:
		figtitle = "OrientationHistogram_Masked_Unweighted_slice_"+str(sliceidx)
		color=Color.red
	
	fig = Plot(figtitle, "Orientation in Degrees", "Distribution of orientation")
	fig.setColor(color);
	fig.add("line",angles, hist)
	fig.setLineWidth(1);
	fig.setLimits(-90, 90, 0, max(hist)*1.1);
	#fig.show()
	imp_fig=fig.getImagePlus()

	return imp_fig



def createMask(imp,roi):
	""" Creates a binary mask (8bit) width width & height like imp.
	value=255 within the roi (or everywhere if roi is None), value=0 elsewhere
	"""

	width=imp.getWidth()
	height=imp.getHeight()

	mask=imp = IJ.createImage("Mask", "8-bit black", width, height, 1);

	maskp = mask.getProcessor();
	maskp.setColor(255);
	if roi is not None:
		maskp.fill(roi);
	else:
		maskp.fill()

	return mask


# ==================================================================================
# general utils
# ==================================================================================

def getFileList(datadir, extension=None, verbose=False):
	"""returns sorted list with all files (full paths) in datadir. Filtered by extension if it is not None."""
	from os import listdir
	from os.path import isfile, join
	
	if extension is None:
		filelist = [os.path.join(datadir,f) for f in sorted(listdir(datadir)) if isfile(join(datadir, f))]
	else:
		filelist = [os.path.join(datadir,f) for f in sorted(listdir(datadir)) if isfile(join(datadir, f)) and f.endswith(extension)]
	
	if verbose:
		IJ.log("Found nr of files: "+str(len(filelist)))
	
	return filelist


def getIntegerSortedFileList(datadir, extension=None, verbose=False):
	"""Similar to getFileList but it does a special sorting (needed for stitcher output filenames) by the increasing integer values within the filename.
	It takes thereby care of missing 0 padding of the filenames.
	["f1_ch0","f2_ch0","f10_ch0"] is returned as ['f1_ch0', 'f2_ch0', 'f10_ch0']
	"""

	from os import listdir
	from os.path import isfile, join

	if extension is None:
		basenames = [str(f) for f in (listdir(datadir)) if isfile(join(datadir, f))]
	else:
		basenames = [str(f) for f in (listdir(datadir)) if isfile(join(datadir, f)) and f.endswith(extension)]

	#print basenames
	basenames.sort(key=lambda f: int(filter(str.isdigit, f)))

	filelist=[os.path.join(datadir,f) for f in basenames]

	if verbose:
		IJ.log("Found nr of files: "+str(len(filelist)))
	
	return filelist


	

# ==================================================================================
# ==================================================================================

main()
