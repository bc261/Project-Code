import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from treedataloader import loader
import pickle
import pprint
import time
import math
import h5py
import os


def massflucuation(tree_data, times, haloid, snap, cutoff):
    """ A function to compute the mass fluctuation for 3 consecutive halos in the main branch of a tree.

    :param tree_data: The tree data dictionary produced by the Merger Graph.
    :param times: A dictionary of the time elapsed at each snapshot. The key corresponds to the snapshot and the
    value corresponds to the elapsed time.
    :param haloid: The halo ID of the current halo being analysed.
    :param snap: The snapshot ID of this halo.
    :param cutoff: The halo mass cutoff in number of particles. Halos under this mass threshold are skipped.

    :return: fluc: The mass fluctuation.
    """

    # Extract the main branch masses
    try:  # If either descendants or progenitors are not present return a NULL value for later removal
        desc_mass = float(tree_data[snap][str(haloid)]['Desc_nPart'][0])
        prog_mass = float(tree_data[snap][str(haloid)]['Prog_nPart'][0])
    except IndexError:
        return -99999
    
    # Extract the present halo's mass
    current_mass = float(tree_data[snap][str(haloid)]['current_halo_nPart'])

    # Only compute the mass fluctuation for 3 consecutive halos above the mass threshold
    if desc_mass < cutoff or prog_mass < cutoff or current_mass < cutoff:
        
        # If any of these masses fall below the threshold return a NULL value
        return -99999

    # Extract the elapsed time of the snapshot
    current_time = times[snap]

    # Compute the progenitor snapshot ID
    if int(snap) > 10:
        prog_snap = '0' + str(int(snap) - 1)
    else:
        prog_snap = '00' + str(int(snap) - 1)

    # Extract the elapsed time of the progenitor snapshot
    prog_time = times[prog_snap]

    # Compute the descendant snapshot ID
    if int(snap) > 8:
        desc_snap = '0' + str(int(snap) + 1)
    else:
        desc_snap = '00' + str(int(snap) + 1)

    # Extract the elapsed time of the descendant snapshot
    desc_time = times[desc_snap]

    # Compute the logarithmic mass growth for each step (prog to present, present to desc)
    logm_desc = (desc_time+current_time)*(desc_mass-current_mass)/((desc_time-current_time)*(desc_mass+current_mass))
    logm_prog = (current_time+prog_time)*(current_mass-prog_mass)/((current_time-prog_time)*(current_mass+prog_mass))

    # Compute the mass fluctuation from the logarithmic mass growths
    fluc = (math.atan(logm_desc) - math.atan(logm_prog))/np.pi

    return fluc

def mainmfluc(tree_data, snap_data, cutoff=10):
    """ A function to compute and plot the mass fluctuation of all sets of 3 halos in a tree'a main branch
    which fall above the mass threshold.
    
    :param tree_data: The tree data dictionary produced by the Merger Graph.
    :param snap_data: The snapshot data dictionary produced by treeloader.py containing all metadata included
    with each snapshot in the simulation.
    :param cutoff: The halo mass cutoff in number of particles. Halos under this mass threshold are skipped.

    :return: None
    """

    # Extract a list of the z=0 halos found in the Merger Graph
    z0_halos = list(tree_data['061'].keys())

    # Remove all halos below the cutoff mass from the list of z=0 halos
    for ind, z0halo in enumerate(z0_halos):
        if tree_data['061'][str(z0halo)]['current_halo_nPart'] < cutoff:
            z0_halos.pop(ind)

    # Create a snapshot list (past to present day) for looping
    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    # Extract snapshot elapsed times and assign to simpler (not nested) dictionary
    times = {}
    for snap in snaplist:
        times[snap] = snap_data[snap]['time']

    # Initialise the list for storing results
    # *** NOTE: This array is initialised with considerably more entries than necessary (namely enough
    # entries for every particle to have a logarithmic mass growth), unused entries are removed after all values
    # have been computed.
    fluc = np.zeros(20000000, dtype=float)

    # Assign the number of halos at z=0 to a variable for progress reporting
    size = len(z0_halos)

    # Initialise the progress counter
    progress = -1

    # Initialise the index counter
    count = -1

    # Loop over all z=0 halos above the mass threshold
    for num, z0halo in enumerate(z0_halos):

        # Print progress if it differs from the previous computed progress
        previous_progress = progress
        progress = int(num/size * 100)
        if progress != previous_progress:
            print('Progress: ', progress , '%')

        # Initialise the number of progenitors
        nprog = -1

        # Initialise the snapshot counter at z=0
        snap = '061'

        # Initialise the main halo pointer
        main = z0halo

        # Walk down the main branch until no progenitors are found (ending the tree)
        while nprog != 0:

            # Increment the index counter
            count += 1

            # Extract the number of progenitors and descendants
            nprog = tree_data[snap][str(main)]['nProg']
            ndesc = tree_data[snap][str(main)]['nDesc']

            # Extract the halo mass to test if the halo falls above the cut off mass
            halo_mass = tree_data[snap][str(main)]['current_halo_nPart']

            # Compute the mass fluctuation if there are descendants, the present halo is above the
            # mass threshold and the current snapshot is not at z=0 since there are no snapshots beyond
            # z=0 meaning no descendants
            if int(snap) < 61 and ndesc != 0 and halo_mass >= cutoff:
                fluc[count] = massflucuation(tree_data, times, main, snap, cutoff=cutoff)

            # If one of the above conditions is false assign a NULL value
            else: count -= 1

            # Find the next halo in the main branch
            if nprog != 0:
                main = tree_data[snap][str(main)]['Prog_haloIDs'][0]

            # Decrement the snapshot counter
            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    # Remove unused array entries
    fluc = fluc[:count+1]

    # Remove instances where progenitors or descendants were missing
    fluc = fluc[np.where(fluc != -99999)]

    # Sort logM results into a histogram
    H, bin_edges = np.histogram(fluc, bins=50)

    # =============== Plot the results ===============

    # Set up figure and axis
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot results
    ax.plot(bin_edges[:-1], H+1, color='r')
    plt.title('Mass fluctuation for Linking length = 0.1')
    
    # Set y axis scale to logarithmic
    ax.set_yscale('log')

    # Label axes
    ax.set_xlabel(r'$\epsilon \pi^{-1}$')
    ax.set_ylabel(r'$N+1$')

    # Save the plot
    plt.savefig('Plots/massFluc_ll=0.1_' + str(cutoff) + '.png', dpi=fig.dpi)

    # Close the plot to avoid memory use issues when looping over many statistics. (In testing it appears
    # matplotlib doesn't completely clear from memory if lots of plots are made in an outer loop)
    plt.close()

    return

def mainmflucComp (tree_data_1, tree_data_2, snap_data, cutoff=10):

    z0_halos_1 = list(tree_data_1['061'].keys())
    z0_halos_2 = list(tree_data_2['061'].keys())
      
    for ind, z0halo in enumerate(z0_halos_1):
        if tree_data_1['061'][str(z0halo)]['current_halo_nPart'] < cutoff:
            z0_halos_1.pop(ind)
    for ind, z0halo in enumerate(z0_halos_2):
        if tree_data_2['061'][str(z0halo)]['current_halo_nPart'] < cutoff:
            z0_halos_2.pop(ind)

    snaplist = []
    for snap in range(0, 62):
        if snap < 10:
            snaplist.append('00' + str(snap))
        elif snap >= 10:
            snaplist.append('0' + str(snap))

    times = {}
    for snap in snaplist:
        times[snap] = snap_data[snap]['time']

    fluc_1 = np.zeros(20000000, dtype=float)
    fluc_2 = np.zeros(20000000, dtype=float)

    size_1 = len(z0_halos_1)
    size_2 = len(z0_halos_2)

    progress = -1

    count = -1

    for num, z0halo in enumerate(z0_halos_1):

        previous_progress = progress
        progress = int(num/size_1 * 100)
        if progress != previous_progress:
            print('Progress: ', progress, '%')

        nprog_1 = -1

        snap = '061'

        main_1 = z0halo
                
        while nprog_1 != 0:

            count += 1
            
            nprog_1 = tree_data_1[snap][str(main_1)]['nProg']
            ndesc_1 = tree_data_1[snap][str(main_1)]['nDesc']
        
            halo_mass_1 = tree_data_1[snap][str(main_1)]['current_halo_nPart']
            
            if int(snap) < 61 and ndesc_1 != 0 and halo_mass_1 >= cutoff:
                fluc_1[count] = massflucuation(tree_data_1, times, main_1, snap, cutoff=cutoff)

            else: count -= 1

            if nprog_1 != 0:
                main_1 = tree_data_1[snap][str(main_1)]['Prog_haloIDs'][0]

            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    fluc_1 = fluc_1[:count+1]
    fluc_1 = fluc_1[np.where(fluc_1 != -99999)]

    count = -1
    
    for num, z0halo in enumerate(z0_halos_2):

        previous_progress = progress
        progress = int(num/size_2 * 100)
        if progress != previous_progress:
            print('Progress: ', progress, '%')

        nprog_2 = -1

        snap = '061'

        main_2 = z0halo
                
        while nprog_2 != 0:

            count += 1
            
            nprog_2 = tree_data_2[snap][str(main_2)]['nProg']
            ndesc_2 = tree_data_2[snap][str(main_2)]['nDesc']
        
            halo_mass_2 = tree_data_2[snap][str(main_2)]['current_halo_nPart']
            
            if int(snap) < 61 and ndesc_2 != 0 and halo_mass_2 >= cutoff:
                fluc_2[count] = massflucuation(tree_data_2, times, main_2, snap, cutoff=cutoff)

            else: count -= 1

            if nprog_2 != 0:
                main_2 = tree_data_2[snap][str(main_2)]['Prog_haloIDs'][0]

            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    fluc_2 = fluc_2[:count+1]
    fluc_2 = fluc_2[np.where(fluc_2 != -99999)]
   
    H_1, bin_edges_1 = np.histogram(fluc_1, bins=50)
    H_2, bin_edges_2 = np.histogram(fluc_2, bins=50)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bin_edges_1[:-1], H_1+1, color='r', label='Linking length = 0.1')
    ax.plot(bin_edges_2[:-1], H_2+1, color='g', label='Linking length = 0.2')
    plt.title('Mass fluctuation Comparison')
    
    ax.set_yscale('log')
    ax.grid(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    ax.set_xlabel(r'$\epsilon \pi^{-1}$')
    ax.set_ylabel(r'$N+1$')

    plt.savefig('Plots/massFlucComp_' + str(cutoff) + '.png', dpi=fig.dpi)

    plt.close()

    return

