"""
"""

import numpy as np
import matplotlib.pyplot as plt


def spot_to_voxel_coords(spot_coords, lb_res):
    """
    voxel size lb_res (x, y, z) e.g. s3 is [0.23*8, 0.23*8, 0.42*4]
    """
    voxel_coords = np.round(spot_coords/lb_res).astype(int)-1
    
    return voxel_coords
    
def filter_raw_spots(spots_rc, intn_th, lb_res, lb):
    """
    lb - image (z, y, x)
    
    apply intensity filter - some thresholds
    apply spatial filter - bounded inside an image
    """
    zlim, ylim, xlim = lb.shape # z, y, x !!!
    
    n_init_spots = len(spots_rc)
    spots_int = spots_rc[:,3]
    spots_xyz = spots_rc[:,:3]
    spots_xyz_v = spot_to_voxel_coords(spots_xyz, lb_res)
    
    # intensity filter
    intn_cond = np.logical_and(
        spots_int > intn_th,
        spots_int < np.percentile(spots_int, 99),
    )

    # spatial filter - bounded within the 3D volume
    spatial_cond = np.all([
        spots_xyz_v[:,0] > 0,
        spots_xyz_v[:,1] > 0,
        spots_xyz_v[:,2] > 0,
        spots_xyz_v[:,0] < xlim,
        spots_xyz_v[:,1] < ylim,
        spots_xyz_v[:,2] < zlim,
    ], axis=0)
    
    filter_cond = np.logical_and(intn_cond, spatial_cond)
    
    print(f"frac. intensity filter: {intn_cond.sum()   /n_init_spots:.2f}")
    print(f"frac. space filter: {spatial_cond.sum()/n_init_spots:.2f}") 
    print(f"frac. combined filter: {filter_cond.sum() /n_init_spots:.2f}")
    
    # remove outside range
    spots_rc = spots_rc[filter_cond]
    
    return spots_rc

def spots_incells_metrics(spots_rc, lb, lb_res):
    """
    spots_rc  - raw spots - filtered 
    !!lb - (z, y, x)  - masks; 0 means outside cells
    lb_res - voxel resolution - (x, y, z)
    """
    spots_xyz = spots_rc[:,:3]
    spots_xyz_v = spot_to_voxel_coords(spots_xyz, lb_res)
    
    spots_incell = (lb[spots_xyz_v[:,2], 
                       spots_xyz_v[:,1], 
                       spots_xyz_v[:,0],] != 0) # z, y, x
    spots_total = len(spots_incell)
    
    space_incell = (lb != 0)
    space_total = lb.size
    
    print(f"frac spots in cells: {np.sum(spots_incell)/spots_total:.2f}")
    print(f"frac voxel in cells: {np.sum(space_incell)/space_total:.2f}")
    
    return 

def scramble_cell_masks(lb):
    """permute within all cell space (non zero)
    """
    # voxels in cell masks
    i, j, k = lb.nonzero()
    # cell labels
    v = lb[i,j,k]
    # shuffled cell labels
    np.random.shuffle(v)
    # shuffled cell masks
    lb_shuff = np.zeros(lb.shape)
    lb_shuff[i,j,k] = v
    
    return lb_shuff

def masks_to_labeled_masks(msk, labeled_cells):
    unq, inv = np.unique(msk.reshape(-1,), return_inverse=True)

    for i in unq:
        if i not in labeled_cells:
            unq[i] = 0

    labeled_msk = unq[inv].reshape(msk.shape)
    
    return labeled_msk

def plot_frac(data, shff, bins=np.arange(0,11,1)):
    """
    """
    fig, ax = plt.subplots()
    
    cnts_data, bins = np.histogram(data, bins)
    cnts_shff, bins = np.histogram(shff, bins)
    
    ax.plot(bins[1:], cnts_shff/cnts_data, '-o')
    ax.set_yscale('log')
    ax.set_yticks([1,0.1,0.05,0.01])
    ax.set_yticklabels([1,0.1,0.05, 0.01])
    
    return fig

def plot_reverse_cumsum(counts, bins=np.arange(0,11,1), ymax=None):
    """
    """
    fig, ax = plt.subplots(figsize=(8,6))
    ax2 = ax.twinx()
    ax.set_xlabel('num spots (normalized)')
    ax.set_ylabel('num cells (cumulative)')
    ax2.set_ylabel('fraction of cells')
    
    n = len(counts)
    cnts, _ = np.histogram(counts, bins)
    rev_cumsum = n-np.cumsum(cnts)
        
    ax.plot(bins[1:], rev_cumsum, '-o', )
    ax2.plot(bins[1:], rev_cumsum/n, '-o', )
    
    if ymax:
        ax.set_ylim(ymin=0, ymax=ymax)
        ax2.set_ylim(ymin=0, ymax=ymax/n)
        
    ax.grid(False)
    ax2.grid(False)
    
    return fig

def plot_reverse_cumsum_complex(bins,
                                counts_list,
                                label_list=None, 
                                color_list=None, 
                                ymax=None):
    """
    """
    fig, axs = plt.subplots(1,2,figsize=(2*6,4))
    ax = axs[0]
    ax2 = ax.twinx()
    ax.set_xlabel('num spots (normalized by mean cell vol)')
    ax.set_ylabel('num cells (cumulative)')
    ax2.set_ylabel('fraction of cells')
    
    if label_list is None:
        label_list = np.arange(len(counts_list))
    if color_list is None:
        color_list = sns.color_palette(n_colors=3)
        
    cumsum_list = []
    n = len(counts_list[0])
    for counts in counts_list:
        assert n == len(counts) # assumes len(counts) is the same
        cnts, _ = np.histogram(counts, bins)
        rev_cumsum = n-np.cumsum(cnts)
        cumsum_list.append(rev_cumsum)
        
    for revcnts, lb, color in zip(cumsum_list, label_list, color_list):
        ax.plot(bins[1:], revcnts, '-', label=lb, color=color, markersize=1)
        ax2.plot(bins[1:], revcnts/n, '-', label=lb, color=color, markersize=1)
    
    if ymax:
        ax.set_ylim(ymin=0, ymax=ymax)
        ax2.set_ylim(ymin=0, ymax=ymax/n)
        
    ax.grid(False)
    ax2.grid(False)
    
    ax = axs[1]
    lb = label_list[1]
    color = color_list[1]
    fdrs = cumsum_list[1]/cumsum_list[0]
    fdrs = np.clip(fdrs, 1e-10, 1,) # keep it in range 1e-10 ~ 1
    ax.plot(bins[1:], -np.log10(fdrs), '-', label=lb, color=color, markersize=1)

    # ax.set_yscale('log')
    # ax.set_yticks([1,0.1,0.05,0.01])
    # ax.set_yticklabels([1,0.1,0.05, 0.01])
    
    ax.set_ylabel('-log10(eFDR) (shuff/data)')
    ax.set_xlabel('num spots (normalized)')
    
    fig.subplots_adjust(wspace=0.4)

    
    return fig, axs, fdrs

