import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
from treedataloader import loader
from smtloader import opentree
import pprint
import time
import h5py
import os
from scipy import stats

def directprogdeschist(tree_data, cutoff):
    """ A function that extracts the number of progenitors and descendants for all halos and
    produces two histograms, one for the number of progenitors and one for the number of
    descendants.

    :param tree_data: The tree data dictionary produced by the Merger Graph.
    :param SMTtreepath: The file path to the SMT algorithms' txt files.

    :return: None
    """

    # Initialise the arrays to store the number of progenitors and descendants
    # *** NOTE: These arrays are initialised with considerably more entries than necessary (namely enough
    # entries for every particle to have a logarithmic mass growth), unused entries are removed after all values
    # have been computed.
    nprogs = np.zeros(20000000, dtype=int)
    ndescs = np.zeros(20000000, dtype=int)

    # Create a snapshot list (past to present day) for looping
    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    # Initialise array index counter
    ind_count = -1
    
    size = len(snaplist)

    progress = -1
    
    # Loop through Merger Graph data assigning each value to the relevant list
    for snap in snaplist:

        # Print the snapshot for progress tracking
        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size * 100)
        if progress != previous_progress:
            print('Progress: ', progress, '%')
  
            
        # Loop over halos within snapshot
        for halo in tree_data[snap].keys():
            
            # Only perform the computation for halos above the mass threshold
            if tree_data[snap][halo]['current_halo_nPart'] >= cutoff:
                continue

            # Increment the index counter
            ind_count += 1
                    
            # Assign the number of progenitors and descendants
            nprogs[ind_count] = tree_data[snap][halo]['nProg']
            ndescs[ind_count] = tree_data[snap][halo]['nDesc']


    # Remove unused entries
    nprogs = nprogs[:ind_count+1]
    ndescs = ndescs[:ind_count+1]

    # Histogram the results with the number of bins defined such that every number between 0 and
    # the max number of progenitors and descendants has a bin.
    prog_bins = int(np.max(nprogs))
    desc_bins = int(np.max(ndescs))
    HprogDMLJ, binedgesprogDMLJ = np.histogram(nprogs, bins=prog_bins)
    HdescDMLJ, binedgesdescDMLJ = np.histogram(ndescs, bins=desc_bins)

    # =============== Plot the results ===============

    # Set up figure and axes
    fig = plt.figure()
    gs = gridspec.GridSpec(2,2)
    gs.update(wspace=0.5, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,:])
    ax2 = fig.add_subplot(gs[1,:])

    # Plot the histograms
    ax1.bar(binedgesprogDMLJ[:-1], HprogDMLJ, color='r', width=0.25, label='Linking_Length=0.1')
    ax2.bar(binedgesdescDMLJ[:-1], HdescDMLJ, color='r', width=0.25, label='Linking_Length=0.1')

    # # Find, histogram and plot the histograms for the SMT softwares
    # # Find the SMT files for each software
    # treefiles = []
    # for root, dirs, files in os.walk(SMTtreepath):
    #     for name in files:
    #         treefiles.append(os.path.join(root, name))
    #
    # # Compute main branch length and histogram for each SMT software
    # for treefile, color in zip(treefiles, plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(list(treefiles))]):
    #
    #     # Get the snapshot progenitor information
    #     try:
    #         software, tot_nprog, prog_num, halo_progs = opentree(treefile)
    #     except UnicodeDecodeError:
    #         continue
    #     print(software)
    #     # Loop through the halo IDs collecting all progenitor numbers
    #     nprogs = np.zeros(len(prog_num.keys()))  # initialise nprog array
    #     for ind, halo in enumerate(prog_num.keys()):
    #
    #         nprogs[ind] = prog_num[halo]
    #
    #     # Histogram the results
    #     prog_bins = int(np.max(nprogs))
    #     Hprog, binedgesprog = np.histogram(nprogs, bins=prog_bins)
    #     ax1.plot(binedgesprog[:-1]+0.5, Hprog, color=color, linestyle='-', label=software)

    # Set y-axis scaling to logarithmic
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.grid(True)
    ax2.grid(True)
    
    # Label axes
    ax1.set_xlabel(r'$N_{Prog}$')
    ax2.set_xlabel(r'$N_{Desc}$')
    ax1.set_ylabel(r'$N$')
    ax2.set_ylabel(r'$N$')
    ax1.set_xlim([0,3])
    ax2.set_xlim([0,3])

    # Include legend
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels)

    # Save the plot as a png
    plt.savefig('Plots/ProgDescNumberHist_ll=0.1_' + str(cutoff) + '_lim_3.png', dpi=fig.dpi)

    return

def progdescsum(tree_data_1, tree_data_2, cutoff):
    nprogs_1 = np.zeros(20000000, dtype=int)
    ndescs_1 = np.zeros(20000000, dtype=int)
    nprogs_2 = np.zeros(20000000, dtype=int)
    ndescs_2 = np.zeros(20000000, dtype=int)

    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    ind_count = -1
    
    size_1 = len(snaplist)

    progress = -1
    
    for snap in snaplist:

        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_1 * 100)
        if progress != previous_progress:
            print('Progress for data set 1: ', progress, '%')

        for halo in tree_data_1[snap].keys():

            if tree_data_1[snap][halo]['current_halo_nPart'] >= cutoff:
                continue

            ind_count += 1

            nprogs_1[ind_count] = tree_data_1[snap][halo]['nProg']
            ndescs_1[ind_count] = tree_data_1[snap][halo]['nDesc']

    nprogs_1 = nprogs_1[:ind_count+1]
    ndescs_1 = ndescs_1[:ind_count+1]

    ind_count = -1
    
    size_2 = len(snaplist)

    progress = -1
    
    for snap in snaplist:

        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_2 * 100)
        if progress != previous_progress:
            print('Progress for data set 2: ', progress, '%')

        for halo in tree_data_2[snap].keys():

            if tree_data_2[snap][halo]['current_halo_nPart'] >= cutoff:
                continue
            
            ind_count += 1

            nprogs_2[ind_count] = tree_data_2[snap][halo]['nProg']
            ndescs_2[ind_count] = tree_data_2[snap][halo]['nDesc']

    nprogs_2 = nprogs_2[:ind_count+1]
    ndescs_2 = ndescs_2[:ind_count+1]

    prog_bins_1 = int(np.max(nprogs_1))
    desc_bins_1 = int(np.max(ndescs_1))
    prog_bins_2 = int(np.max(nprogs_2))
    desc_bins_2 = int(np.max(ndescs_2))

    HprogDMLJ_1, binedgesprogDMLJ_1 = np.histogram(nprogs_1, bins=prog_bins_1)
    HdescDMLJ_1, binedgesdescDMLJ_1 = np.histogram(ndescs_1, bins=desc_bins_1)
    HprogDMLJ_2, binedgesprogDMLJ_2 = np.histogram(nprogs_2, bins=prog_bins_2)
    HdescDMLJ_2, binedgesdescDMLJ_2 = np.histogram(ndescs_2, bins=desc_bins_2)

    HprogDMLJ_01=HprogDMLJ_1.astype(float)
    HdescDMLJ_01=HdescDMLJ_1.astype(float)
    HprogDMLJ_02=HprogDMLJ_2.astype(float)
    HdescDMLJ_02=HdescDMLJ_2.astype(float)

    sumprog_1 = np.cumsum(HprogDMLJ_1[::-1], dtype=np.float64)[::-1]
    sumdesc_1 = np.cumsum(HdescDMLJ_1[::-1], dtype=np.float64)[::-1]
    sumprog_2 = np.cumsum(HprogDMLJ_2[::-1], dtype=np.float64)[::-1]
    sumdesc_2 = np.cumsum(HdescDMLJ_2[::-1], dtype=np.float64)[::-1]

    maxsumprog_1=max(sumprog_1)
    maxsumdesc_1=max(sumdesc_1)
    maxsumprog_2=max(sumprog_2)
    maxsumdesc_2=max(sumprog_2)

    divprog_1 = sumprog_1/maxsumprog_1
    divdesc_1 = sumdesc_1/maxsumdesc_1
    divprog_2 = sumprog_2/maxsumprog_2
    divdesc_2 = sumdesc_2/maxsumdesc_2
    
    fig = plt.figure()
    gs = gridspec.GridSpec(2,2)
    gs.update(wspace=0.5, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,:])
    ax2 = fig.add_subplot(gs[1,:])

    ax1.plot(binedgesprogDMLJ_1[:-1], divprog_1 , color='r', label='Linking Length = 0.1')
    ax1.plot(binedgesprogDMLJ_2[:-1], divprog_2 , color='g', label='Linking Length = 0.2')
    ax2.plot(binedgesdescDMLJ_1[:-1], divdesc_1 , color='r', label='Linking Length = 0.1')
    ax2.plot(binedgesdescDMLJ_2[:-1], divdesc_2 , color='g', label='Linking Length = 0.2')
  
    ax1.grid(True)
    ax2.grid(True)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    #ax1.set_xlim(0,10)
    #ax2.set_xlim(0,10)
    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels)
    
    ax1.set_xlabel(r'Number')
    ax2.set_xlabel(r'Number')
    ax1.set_ylabel(r'Histogram of Progenitors')
    ax2.set_ylabel(r'Histogram of Descendants')

    plt.savefig('Plots/ProgDescHist_' + str(cutoff) + '.png', dpi=fig.dpi)

    return

def progdescstat(tree_data_1, tree_data_2, cutoff):
    nprogs_1 = np.zeros(20000000, dtype=int)
    ndescs_1 = np.zeros(20000000, dtype=int)
    nprogs_2 = np.zeros(20000000, dtype=int)
    ndescs_2 = np.zeros(20000000, dtype=int)

    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    ind_count = -1
    
    size_1 = len(snaplist)

    progress = -1
    
    for snap in snaplist:

        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_1 * 100)
        if progress != previous_progress:
            print('Progress for data set 1: ', progress, '%')

        for halo in tree_data_1[snap].keys():

            if tree_data_1[snap][halo]['current_halo_nPart'] >= cutoff:
                continue
            
            ind_count += 1

            nprogs_1[ind_count] = tree_data_1[snap][halo]['nProg']
            ndescs_1[ind_count] = tree_data_1[snap][halo]['nDesc']
    
    nprogs_1 = nprogs_1[:ind_count+1]
    ndescs_1 = ndescs_1[:ind_count+1]
    
    ind_count = -1
        
    size_2 = len(snaplist)

    progress = -1
    
    # Loop through Merger Graph data assigning each value to the relevant list
    for snap in snaplist:

        # Print the snapshot for progress tracking
        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_2 * 100)
        if progress != previous_progress:
            print('Progress for data set 2: ', progress, '%')

        for halo in tree_data_2[snap].keys():

            if tree_data[snap][halo]['current_halo_nPart'] >= cutoff:
                continue

            ind_count += 1

            nprogs_2[ind_count] = tree_data_2[snap][halo]['nProg']
            ndescs_2[ind_count] = tree_data_2[snap][halo]['nDesc']

    nprogs_2 = nprogs_2[:ind_count+1]
    ndescs_2 = ndescs_2[:ind_count+1]

    KS_test_prog = stats.ks_2samp(nprogs_1, nprogs_2)
    KS_test_desc = stats.ks_2samp(ndescs_1, ndescs_2)
    
    print('KS_test_prog', KS_test_prog)
    print('KS_test_desc', KS_test_desc)

    return

def progdeschistComp(tree_data_1, tree_data_2, cutoff):
    nprogs_1 = np.zeros(20000000, dtype=int)
    ndescs_1 = np.zeros(20000000, dtype=int)
    nprogs_2 = np.zeros(20000000, dtype=int)
    ndescs_2 = np.zeros(20000000, dtype=int)

    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    ind_count = -1
    
    size_1 = len(snaplist)

    progress = -1
    
    for snap in snaplist:

        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_1 * 100)
        if progress != previous_progress:
            print('Progress for data set 1: ', progress, '%')
  
        for halo in tree_data_1[snap].keys():
            
            if tree_data_1[snap][halo]['current_halo_nPart'] >= cutoff:
                continue

            ind_count += 1
            
            nprogs_1[ind_count] = tree_data_1[snap][halo]['nProg']
            ndescs_1[ind_count] = tree_data_1[snap][halo]['nDesc']

    nprogs_1 = nprogs_1[:ind_count+1]
    ndescs_1 = ndescs_1[:ind_count+1]

    ind_count = -1

    size_2 = len(snaplist)

    progress = -1
    
    for snap in snaplist:

        int_snap = int(snap)
        previous_progress = progress
        progress = int(int_snap/size_2 * 100)
        if progress != previous_progress:
            print('Progress for data set 2: ', progress, '%')
  
        for halo in tree_data_2[snap].keys():
            
            if tree_data_2[snap][halo]['current_halo_nPart'] >= cutoff:
                continue

            ind_count += 1
            
            nprogs_2[ind_count] = tree_data_2[snap][halo]['nProg']
            ndescs_2[ind_count] = tree_data_2[snap][halo]['nDesc']

    nprogs_2 = nprogs_2[:ind_count+1]
    ndescs_2 = ndescs_2[:ind_count+1]
    
    prog_bins_1 = int(np.max(nprogs_1))
    desc_bins_1 = int(np.max(ndescs_1))
    prog_bins_2 = int(np.max(nprogs_2))
    desc_bins_2 = int(np.max(ndescs_2))

    HprogDMLJ_1, binedgesprogDMLJ_1 = np.histogram(nprogs_1, bins=prog_bins_1)
    HdescDMLJ_1, binedgesdescDMLJ_1 = np.histogram(ndescs_1, bins=desc_bins_1)
    HprogDMLJ_2, binedgesprogDMLJ_2 = np.histogram(nprogs_2, bins=prog_bins_2)
    HdescDMLJ_2, binedgesdescDMLJ_2 = np.histogram(ndescs_2, bins=desc_bins_2)

    fig = plt.figure()
    gs = gridspec.GridSpec(2,2)
    gs.update(wspace=0.5, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,:])
    ax2 = fig.add_subplot(gs[1,:])

    ax1.bar(binedgesprogDMLJ_1[:-1], HprogDMLJ_1 , color='r', label='Linking Length = 0.1')
    ax2.bar(binedgesdescDMLJ_1[:-1], HdescDMLJ_1 , color='r', label='Linking Length = 0.1')
    ax1.bar(binedgesprogDMLJ_2[:-1], HprogDMLJ_2 , color='b', align='edge', width=0.4, label='Linking Length = 0.2')
    ax2.bar(binedgesdescDMLJ_2[:-1], HdescDMLJ_2 , color='b', align='edge', width=0.4, label='Linking Length = 0.2')

    ax1.grid(True)
    ax2.grid(True)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
    ax1.set_xlabel(r'$N_{Prog}$')
    ax2.set_xlabel(r'$N_{Desc}$')
    ax1.set_ylabel(r'$N$')
    ax2.set_ylabel(r'$N$')
    #ax1.set_xlim([0,3])
    #ax2.set_xlim([0,3])

    # Include legend
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels)

    # Save the plot as a png
    plt.savefig('Plots/ProgDescHistComp_' + str(cutoff) + '.png', dpi=fig.dpi)

    return
