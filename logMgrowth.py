import numpy as np
import matplotlib.pyplot as plt


def logMassGrowth(tree_data, times, haloid, snap, cutoff):
    """ A function that computes the slope of the mass history (logarithmic mass growth) for a
    single halo and it's main descendant.

    :param tree_data: The tree data dictionary produced by the Merger Graph.
    :param times: A dictionary of the time elapsed at each snapshot. The key corresponds to the snapshot and the
    value corresponds to the elapsed time.
    :param haloid: The halo ID of the current halo being analysed.
    :param snap: The snapshot ID of this halo.
    :param cutoff: The halo mass cutoff in number of particles. Halos under this mass threshold are skipped.

    :return: logm: The logarithmic mass growth.
    """
    # Extract the mass of the halo and the main descendant (number of particles contained in each) from the
    # Merger Graph dictionary
    desc_mass = tree_data[snap][str(haloid)]['Desc_nPart'][0]
    current_mass = tree_data[snap][str(haloid)]['current_halo_nPart']

    # Only compute the logarithmic mass growth for 2 consecutive halos above the mass threshold
    if desc_mass < cutoff or current_mass < cutoff:
        return -99999

    # Extract the elapsed time of the snapshot
    current_time = times[snap]

    # Compute the descendant snapshot ID
    if int(snap) > 8:
        desc_snap = '0' + str(int(snap) + 1)
    else:
        desc_snap = '00' + str(int(snap) + 1)

    # Extract the descendant snapshot's elapsed time
    desc_time = times[desc_snap]

    # Compute the logarithmic mass growth
    logm = (float((desc_time+current_time)*(desc_mass-current_mass)) /
            float((desc_time-current_time)*(desc_mass+current_mass)))

    return logm


def mainLogM(tree_data, snap_data, cutoff=10):
    """ A function which computes the slope of the mass history for all pairs of halo and
    descendant halo found in the Merger Graph and produces a histogram.

    :param tree_data: The tree data dictionary produced by the Merger Graph.
    :param snap_data: The snapshot data dictionary produced by treeloader.py containing all metadata included
    with each snapshot in the simulation.
    :param cutoff: The halo mass cutoff in number of particles. Halos under this mass threshold are skipped.

    :return: None
    """

    # Extract a list of the halos present at z=0
    z0_halos = list(tree_data['061'].keys())

    # Remove all halos below the cutoff mass by removing their keys from the z=0 halos list
    for ind, z0halo in enumerate(z0_halos):
        if tree_data['061'][str(z0halo)]['current_halo_nPart'] < cutoff:
            z0_halos.pop(ind)

    # Create the snapshot list (past to present day) for looping
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
    logM = np.zeros(20000000, dtype=float)

    # Assign the number of halos at z=0 to a variable for progress reporting
    size = len(z0_halos)

    # Initialise the progress counter
    progress = -1

    # Initialise an index counter to keep track of how many logarithmic mass growth array entires have been filled
    ind_count = -1

    # Loop over the z=0 halos
    for num, z0halo in enumerate(z0_halos):

        # Compute and print the progress
        previous_progress = progress
        progress = int(num/size * 100)
        if progress != previous_progress:
            print('Progress: ', progress, '%')

        # Initialise the number of progenitors for the while loop condition
        nprog = -1

        # Initialise the snapshot counter at z=0
        snap = '061'

        # Initialise the main halo pointer
        main = z0halo

        # Loop until the main branch of the current z=0 halo is completed (no progenitors are found in the
        # next snapshot)
        while nprog != 0:

            # Extract the number of progentiors and descendants of the current halo
            nprog = tree_data[snap][str(main)]['nProg']
            ndesc = tree_data[snap][str(main)]['nDesc']

            # Extract halo mass to test if the halo falls above the cut off mass
            halo_mass = tree_data[snap][str(main)]['current_halo_nPart']

            # Avoid computing logM on the first step since there are no descendants
            if int(snap) < 61 and ndesc != 0 and halo_mass >= cutoff:

                # Increment index counter
                ind_count += 1

                # Compute the logarithmic mass growth
                logM[ind_count] = logMassGrowth(tree_data, times, main, snap, cutoff)

            # Find the next halo in the main branch (progenitor IDs are mass ordered)
            if nprog != 0:
                main = tree_data[snap][str(main)]['Prog_haloIDs'][0]

            # Decrement the snapshot counter
            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    # Remove unused entries in the array using the used index counter
    logM = logM[:ind_count+1]

    # Remove the filler value of -99999 used when halos aren't consecutively above the mass threshold
    logM = logM[np.where(logM != -99999)]

    # Sort logM results into a histogram with high resolution bins
    H, bin_edges = np.histogram(logM, bins=1000)

    # =============== Plot the results ===============

    # Set up the figure and axis
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot the results with arctan scaling of the x axis
    ax.plot(np.arctan(bin_edges[:-1]), H+1, color='r')

    # Set y scale to be logarithmic
    ax.set_yscale('log')
    ax.grid(True)
    
    # Label axes
    ax.set_xlabel(r'$\alpha_{M}$')
    ax.set_ylabel(r'$N+1$')

    # Save the plot
    plt.savefig('Plots/logMgrowth_ll=0.1_' + str(cutoff) + '.png', dpi=fig.dpi)

    # Close the plot to avoid memory use issues when looping over many statistics. (In testing it appears
    # matplotlib doesn't completely clear from memory if lots of plots are made in an outer loop)
    plt.clf()

    return

def mainLogMComp(tree_data_1, tree_data_2, snap_data, cutoff=10):

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

    logM_1 = np.zeros(20000000, dtype=float)
    logM_2 = np.zeros(20000000, dtype=float)

    size_1 = len(z0_halos_1)
    size_2 = len(z0_halos_2)

    progress = -1

    ind_count = -1

    for num, z0halo in enumerate(z0_halos_1):

        previous_progress = progress
        progress = int(num/size_1 * 100)
        if progress != previous_progress:
            print('Progress for data set 1: ', progress, '%')

        nprog_1 = -1

        snap = '061'

        main_1 = z0halo
                
        while nprog_1 != 0:

            nprog_1 = tree_data_1[snap][str(main_1)]['nProg']
            ndesc_1 = tree_data_1[snap][str(main_1)]['nDesc']
        
            halo_mass_1 = tree_data_1[snap][str(main_1)]['current_halo_nPart']
            
            if int(snap) < 61 and ndesc_1 != 0 and halo_mass_1 >= cutoff:

                ind_count += 1

                logM_1[ind_count] = logMassGrowth(tree_data_1, times, main_1, snap, cutoff)
               
            if nprog_1 != 0:
                main_1 = tree_data_1[snap][str(main_1)]['Prog_haloIDs'][0]

            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    logM_1 = logM_1[:ind_count+1]            
    logM_1 = logM_1[np.where(logM_1 != -99999)]      

    ind_count = -1

    progress = -1
    
    for num, z0halo in enumerate(z0_halos_2):

        previous_progress = progress
        progress = int(num/size_2 * 100)
        if progress != previous_progress:
            print('Progress for data set 2: ', progress, '%')

        nprog_2 = -1

        snap = '061'

        main_2 = z0halo
                
        while nprog_2 != 0:

            nprog_2 = tree_data_2[snap][str(main_2)]['nProg']
            ndesc_2 = tree_data_2[snap][str(main_2)]['nDesc']

            halo_mass_2 = tree_data_2[snap][str(main_2)]['current_halo_nPart']
            
            if int(snap) < 61 and ndesc_2 != 0 and halo_mass_2 >=cutoff :

                ind_count += 1

                logM_2[ind_count] = logMassGrowth(tree_data_2, times, main_2, snap, cutoff)

            if  nprog_2 !=0:
                
                main_2 = tree_data_2[snap][str(main_2)]['Prog_haloIDs'][0]

            if int(snap) > 10:
                snap = '0' + str(int(snap) - 1)
            else:
                snap = '00' + str(int(snap) - 1)

    logM_2 = logM_2[:ind_count+1]
    logM_2 = logM_2[np.where(logM_2 != -99999)]

    H_1, bin_edges_1 = np.histogram(logM_1, bins=1000)
    H_2, bin_edges_2 = np.histogram(logM_2, bins=1000)
      
    fig = plt.figure()
    ax = fig.add_subplot(111)
 
    ax.plot(np.arctan(bin_edges_1[:-1]), H_1+1, color='r', label='Linking Length = 0.1')
    ax.plot(np.arctan(bin_edges_2[:-1]), H_2+1, color='g', label='Linking Length = 0.2')
    
    ax.set_yscale('log')
    ax.grid(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    ax.set_xlabel(r'$\alpha_{M}$')
    ax.set_ylabel(r'$N+1$')

    plt.savefig('Plots/logMgrowthComp_' + str(cutoff) + '.png', dpi=fig.dpi)

    plt.clf()

    return
