# All functions used in the jupyter notebook, to avoid cluttering the notebook.

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_hist(data,slice_selection, use_weighted, normalize_by_maxval=False,save_dir=None, color_map=None, ylim=None):    
    """Plots orientation distribution in different flavors.
    data: numpy array with (rows,cols)=(slices,histogram count (180bins)) (obtained from fiji script + processed in notebook)
    slice_selection: None or a list of indices, e.g. [1] or [2,5]. one-based
    use_weighted. boolean. Whether a coherency-weighted histogram is used. Must match the data! Here used for correct
            axis labeling.
    normalize_by_maxval: boolean. whether to normalize the distribution of each slice individually by its max,  such 
            that the max value (most frequent hist count) is =1. 
    save_dir: if not None: saves plots into this directory. Folder must exist on disk.
    color_map: None (uses "copper") or any of these strings: https://matplotlib.org/tutorials/colors/colormaps.html
    ylim: if not None: sets the ylim
    """
    
    # prepare save name for figure
    figname="OrientationHistogram"
    descr1="_AllSlices" if (slice_selection is None) else "_SelectedSlices"
    descr2="_CoherencyWeightedCount" if (use_weighted) else "_UnweightedCount"
    descr3="_NormalizedToMax1" if (normalize_by_maxval) else ""
    figname=figname+descr1+descr2+descr3+".png"
    
    # prepare figure description
    if use_weighted:
        if normalize_by_maxval:
            ylabel="frequency (weighted, normalized by max)"
            title="Orientation distribution (weighted with coherencies, per slice normalized by max)"
        else:
            ylabel="frequency (weighted)"
            title="Orientation distribution (weighted with coherencies)"
    else:
        if normalize_by_maxval:
            ylabel="frequency (normalized by max)"
            title="Orientation distribution (unweighted, per slice normalized by max)"
        else:
            ylabel="frequency"
            title="Orientation distribution (unweighted)"
      
    # init plot
    angles=np.linspace(-89.5,89.5,180) # center of bins
    plt.figure(figsize=(12,7))
    ax = plt.subplot(111)
    
    # use all slices if not specified
    if slice_selection is None:
        slice_selection=np.arange(1,1+data.shape[0])
    
    # prepare array with colors
    if color_map is None:
        color_map="copper"
    colors=np.squeeze(plt.get_cmap(color_map)([np.linspace(0,1,len(slice_selection))]))
    
    # if only one slice plotted, need to make colors array 2d
    if colors.ndim==1:
        colors=np.expand_dims(colors,axis=0)
        
    # actual plotting
    for idx in range(len(slice_selection)):
        s=slice_selection[idx]
        y=data[s-1,:] # slices are one-based

        if normalize_by_maxval:
            maxval=np.max(y)
            y=y/maxval

        plt.plot(angles,y,color=colors[idx],label="slice "+str(s)) 
        
    plt.xlim([-90,90])
    plt.xlabel("orientation (in degree)",size=14)
    plt.ylabel(ylabel,size=14)
    plt.title(title,size=14)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    # y-axis adjustement
    if ylim is not None:
        plt.ylim(ylim)
    
    # save the figure
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir,figname))

    plt.show()
    
   
    
def plot_coherency(data,use_weighted, save_dir=None):       
    if not use_weighted:
        print("Cannot calculate average coherency with unweighted histogram. Rerun the notebook with 'use_weighted=True'")
        return
        
    # unweighted histogram sums to =1. therefore the weighted sums to mean(coherency)
    avg_coherencies=np.sum(data,1) 
    xax=np.arange(1,len(avg_coherencies)+1)

    plt.figure(figsize=(9,5))
    plt.plot(avg_coherencies,".-")
    plt.ylim([0,1])
    plt.xlabel("slice id",size=16)
    plt.ylabel("average coherency",size=16)
    plt.title("average coherency (measure for orientation confidence) per slice",size=16)
    plt.xticks(xax)
    
    # save the figure
    if save_dir is not None:
        fn=os.path.join(save_dir,"Average_Coherency_per_Slice.png")
        plt.savefig(fn)
        
    plt.show()

  
    
def get_filename_csv_histdata(data_folder, use_weighted):
    """Returns the full path filename of the csv file which contains the histogram data, takes thereby care of 
    whether the weighted/unweighted data is being used.
    Basenames of the csv files are hardcoded.
    """
    fn_weightedhist="orientationDistribution_CoherencyWeighted.csv"
    fn_unweightedhist="orientationDistribution_Unweighted.csv"

    if use_weighted:
        fn_load=os.path.join(data_folder,fn_weightedhist)
    else:
        fn_load=os.path.join(data_folder,fn_unweightedhist)

    if not os.path.isfile(fn_load):
        print("ERROR: Cannot find the csv file with histogram data: \n",fn_load)
    else:
        print("Found csv file with histogram data:",fn_load)

    return fn_load
    
    
def quick_plot_of_peak_selection(datasets,use_weighted,sliceids_longitudinal, sliceids_circular, titles, save_dir=None, ylim=None):    
    """ Quick visualization to check whether the selected circular/longitudinal slices (peaks) are correct.
    Comparable to plot_hist but less fancy (however more convenient colormap for this purpose). 
    Does allow for plotting multiple datasets at once.
    data: list[arrays], or a single numpy array with (rows,cols)=(slices,histogram count (180bins)) (obtained from fiji script + processed in notebook)
    use_weighted: for axis labeling
    sliceids_longitudinal: list[ints] or single int of the slice id with the longitudinal fiber orientation peak
    sliceids_circular: see above
    titles: list[strings] or string with typically the dataset name. 
    save_dir: if not None: saves plots into this directory. Folder must exist on disk.
    ylim: if not None: sets the ylim
    """
    
    # save name for figure
    figname="OverviewPeakSelection.png"
    
    # prepare figure description
    if use_weighted:
        ylabel="frequency (weighted)"
    else:
        ylabel="frequency"
        
    # convert single dataset also to list
    if not isinstance(datasets,list):
        datasets=[datasets]
        sliceids_longitudinal=[sliceids_longitudinal]
        sliceids_circular=[sliceids_circular]
        titles=[titles]
        
    assert len(datasets) == len(sliceids_longitudinal) == len(sliceids_circular) == len(titles)
    
    # bin centers
    angles=np.linspace(-89.5,89.5,180) 
    
    # subplots layout
    ncols=1 if len(datasets)==1 else 2
    nrows=int(1+len(datasets)/ncols)
    
    plt.figure(figsize=(6*ncols,4*nrows))
    plt.subplots_adjust(hspace=0.4) # default=0.2
    
    # one subplot per dataset
    for data_id in range(len(datasets)):
        plt.subplot(nrows,ncols,data_id+1)
            
        data=datasets[data_id]
        slice_long=sliceids_longitudinal[data_id]
        slice_circ=sliceids_circular[data_id]
        title=titles[data_id]
        
        for idx in range(data.shape[0]):
            plt.plot(angles,data[idx,:],color=[0.7,0.7,0.7])
        
        plt.plot(angles,data[slice_long-1,:],color=[0.9,0,0]) #red
        plt.plot(angles,data[slice_circ-1,:],"blue") 
        
        plt.xlim([-90,90])
        plt.xlabel("orientation (in degree)",size=12)
        plt.ylabel(ylabel,size=12)
        plt.title(title+" (slices: "+str(slice_circ)+", "+str(slice_long)+")",size=14)

        # y-axis adjustement
        if ylim is not None:
            plt.ylim(ylim)


    # save the figure
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir,figname))

    plt.show()
    
    


    
def _sliceid_of_highest_peak_within_tolerance(maxvals, dist, tolerance):
    """returns the id of the slice with the highest peak value (maxvals), for which the peak location distance (dist)
    to the expected  max location is less than peak tolerance (in deg)
    slice_id is one based!.
    """
    sub_sliceidxs=np.where(dist<tolerance)[0] # zero-based
    if len(sub_sliceidxs)==0:
        print("ERROR: No peak with suitable location found. Check expected peak locations and peak tolerance.")
        return None
    
    sub_maxvals=maxvals[sub_sliceidxs]
    max_id=sub_sliceidxs[np.argmax(sub_maxvals)]+1 # convert to one-based

    return max_id


def _sliceids_of_longitudinal_and_circular_hist_peak(data, peak_circular_expected, peak_longitudinal_expected, peak_tolerance=20):
    """Finds the slices (one-based) for which the strongest circular, resp. longitudinal orientation exists (i.e.
    the peak in the histogram is highest)
    data: numpy array with (rows,cols)=(nslices,180)
    peak_circular_expected: orientation in degree where the peak of circular fiber layers is expected, typically 0.
                    between [-90,90]
    peak_longitudinal_expected: see above. typically 90. between [-90,90]
    peak_tolerance: all hist peaks within this tolerance (in degree) are considered a longitudinal resp circular layer
    returns: circ_id, long_id: ids (one-based) of the circular and longitudinal layer with highest peak
    """
    # histogram peak per slice
    max_idxs=np.argmax(data,1) # location (as array idx)
    max_values=data[np.arange(data.shape[0]),max_idxs] # peak height

    # location of hist peak in degree
    angles=np.linspace(-89.5,89.5,180) # assumes 180bins
    max_loc_degrees=angles[max_idxs]
    
    # find the slices with the strongest peak for longitudinal & circular fibers   
    def compute_angle_dist(angle, angle_expected):
        """takes care of 180deg ambiguity and maps angle difference into [0,90]"""
        dist=(abs(angle-angle_expected))%180 # in case an angle was >90 or <-90
        dist[dist>90]=abs(180-dist[dist>90])
        return dist # in degree
        
    dist_longitudinal=compute_angle_dist(max_loc_degrees,peak_longitudinal_expected) 
    dist_circular=compute_angle_dist(max_loc_degrees,peak_circular_expected)

    # slice ids of circular/longitudinal peaks
    circ_id=_sliceid_of_highest_peak_within_tolerance(max_values,dist_circular,peak_tolerance)
    long_id=_sliceid_of_highest_peak_within_tolerance(max_values,dist_longitudinal,peak_tolerance)  
    
    return circ_id, long_id



def load_datasets_and_find_peaks(list_of_folders,use_weighted, peak_circular_expected, peak_longitudinal_expected):
    """Processes all folders in list_of_folders: loads for each folder the histogram csv file, then finds the 
    slice-ids (one-based) with the highest peak (for circular and longitudinal fibers separately).
    Params:
    use_weighted: should always be true
    peak_circular_expected, peak_longitudinal_expected: see _sliceids_of_longitudinal_and_circular_hist_peak
    Returns:
    data_list: list of numpy arrays of size (nslices, 180), one entry per folder of list_of_folders
    id_list_circular: list of slice id's (one based) for which the histogram of the circular fibers was highest.
                    One id per dataset.
    id_list_longitudinal: see above   
    """

    data_list=[]
    id_list_circular=[]
    id_list_longitudinal=[]
    
    for datadir in list_of_folders:
        # load data and convert to numpy array
        fn_load=get_filename_csv_histdata(datadir,use_weighted)
        df=pd.read_csv(fn_load)
        data=(df.iloc[:,1:]).values # first/second,.. row: slice 1,2,..

        # find slices with highest peaks
        circ_id, long_id = _sliceids_of_longitudinal_and_circular_hist_peak(data,peak_circular_expected,peak_longitudinal_expected)
        
        data_list.append(data)
        id_list_circular.append(circ_id)
        id_list_longitudinal.append(long_id)
        
    return data_list,id_list_circular,id_list_longitudinal


def get_statistics(datasets, sliceids_long, sliceids_circ):
    """Computes the mean and stddev between the orientation distributions of the different experiemnts.
    For each experiemnt the longitudinal and circular distribution with the highest peak is chosen 
    params: see output of load_datasets_and_find_peaks()
    returns:
    ymean_circ, ystddev_circ, ymean_long, ystddev_long: the mean and stddev of the distributions, array of shape (180,)
    """
    # collect circular, resp. longitudinal max-peak curves 
    y_circ=collect_peak_distributions(datasets,sliceids_circ)
    y_long=collect_peak_distributions(datasets,sliceids_long)

    # get statistics
    ymean_circ=np.mean(y_circ,0)
    ystddev_circ=np.std(y_circ,0)
    ymean_long=np.mean(y_long,0)
    ystddev_long=np.std(y_long,0)    
    
    return ymean_circ, ystddev_circ, ymean_long, ystddev_long

    
def collect_peak_distributions(datasets, sliceids):
    """collect curves with highest peak. (sliceid: one-based). returns array (nsamples,180)"""
    y=np.zeros((len(datasets),datasets[0].shape[-1]))
    for i in range(y.shape[0]):
        y[i]=datasets[i][sliceids[i]-1,:] # one-based
    return y


def plot_average_histogram(ym_circ,ystd_circ,ym_long,ystd_long, numsamples=None, save_dir=None,ylim=None):    
    """ Plots the average +- stddev orientaiton distribution of the orientation histogram (for longitudinal and circular).
    Purpose: plot result of get_statistics().
    Params:
    ym_circ, ystd_circ: mean and stddev of circular fiber distribution
    ym_long, ystd_long: ... longitudinal ...
    [removed: use_weighted]
    numsamples: optional: provide the number of datasets that was used (int) -> added to title
    save_dir: if not None: saves plots into this directory. Folder must exist on disk.
    ylim: if not None: sets the ylim
    """
    
    figname="AveragedOrientationHistogram.png"
    ylabel="frequency (mean+-stddev)"
    title="Average orientation distribution"
    if numsamples is not None:
        title+=" (n="+str(numsamples)+")"
           
    # bin centers
    angles=np.linspace(-89.5,89.5,180) 
       
    plt.figure(figsize=(7,5))
    ax = plt.subplot(111)
    
    # plotting
    plt.plot(angles,ym_circ, 'b-',label="circular")
    plt.fill_between(angles,ym_circ-ystd_circ,ym_circ+ystd_circ,color=[0.8,0.8,0.8]) # light gray

    plt.plot(angles,ym_long, '-',color=[.9,0,0],label="longitudinal") # red
    plt.fill_between(angles,ym_long-ystd_long,ym_long+ystd_long,color=[0.8,0.8,0.8]) 

    plt.xlim([-90,90])
    plt.xlabel("orientation (in degree)",size=14)
    plt.ylabel(ylabel,size=14)
    plt.title(title,size=14)
   
    ax.legend()

    # y-axis adjustement
    if ylim is not None:
        plt.ylim(ylim)

    # save the figure
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir,figname))

    plt.show()