import numpy as np
import matplotlib.pyplot as plt
import h5py
                    
def test(tree_data_1, tree_data_2, snap_data):
    z0_halos_1 = list(tree_data_1['061'].keys())
    z0_halos_2 = list(tree_data_2['061'].keys())

    z0_halos_1 = [int(x) for x in z0_halos_1]
    z0_halos_2 = [int(x) for x in z0_halos_2]
     
    z0_halos_1.sort()
    z0_halos_2.sort()

    print('z0_halos_1:', z0_halos_1)
    print('z0_halos_2:', z0_halos_2)

    return

def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)
    return

def SumMassCompPlot(halo_snapshot):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    #print('f.keys():',list(f.keys()))
    datasets = [item for item in f.values() if isinstance(item, h5py.Dataset)]
    groups = [item for item in f.values() if isinstance(item, h5py.Group)]
    print('len(groups):',len(groups))
    print('len(datasets):',len(datasets))
    base_list=len(f['/Halo_IDs/'])
    nHalo = len(groups)
    print('nHalo:',nHalo)
    main_mass = np.empty(nHalo,dtype=np.int32)
    nHaloSub=np.empty(nHalo,dtype=np.int32)
    progress = -1
    
    # First loop to count halos
    iSub=0
    nSub=0
    for i in range(0,nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress main_IDs: ', progress, '%')
        main_mass[i]=f['/'+str(i)].attrs['halo_nPart']
        nHaloSub[i] = len([item for item in f['/'+str(i)].values() if isinstance(item, h5py.Group)])
    nSub=np.sum(nHaloSub)


    # Array to hold subhalo masses
    sub_mass=np.empty(nSub,dtype=np.int32)
    split_sub_mass = np.empty(nSub,dtype=np.int32)
    tot_sub_mass=np.empty(nHalo,dtype=np.int32)
    iSub=0
    
    progress = -1
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs: ', progress, '%')
        for j in range(0,nHaloSub[i]):
            sub_mass[iSub]=f['/'+str(i)+'/'+str(j)].attrs['subhalo_nPart']
            #print('len(sub_mass):',len(sub_mass))
            #print('sub_mass:',sub_mass)
            iSub +=1

    print('nhaloSub:',nHaloSub)
    print('len(nHaloSub):',len(nHaloSub))
    split_sub_mass = np.split(sub_mass, np.cumsum(nHaloSub)[:-1])
    tot_sub_mass = [sum(x) for x in split_sub_mass]
    #print('sub_mass:',sub_mass)
    #print('len(nHaloSub):',len(main_mass))
    #print('main_mass:',main_mass)
    print('len(main_mass):',len(main_mass))
    print('main_mass:',main_mass)
    print('len(sub_mass):',len(sub_mass))
    print('sub_mass:',sub_mass)
    print('len(split_sub_mass):',len(split_sub_mass))
    #print('split_sub_mass:',split_sub_mass)
    #print('tot_sub_mass:',list(tot_sub_mass))
    print('len(tot_sub_mass):',len(tot_sub_mass))

    fig = plt.figure()
    ax = fig.add_subplot(111)
 
    ax.plot(main_mass, tot_sub_mass,'.', color='r')

    ax.grid(True)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_title('The Sum of 0.1 Linking Lengths')
    ax.set_xlabel(r'$M_{0.2}$')
    ax.set_ylabel(r'$M_{0.1}$')

    plt.savefig('Plots/SumMassComp.png', dpi=fig.dpi)

    plt.clf()
    
    return

def NumOcup(halo_snapshot):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    #print('f.keys():',list(f.keys()))
    groups = [item for item in f.values() if isinstance(item, h5py.Group)]
    print('len(groups):',len(groups))
    nHalo = len(groups)
    nsub = np.empty(nHalo,dtype=np.int32)
    print('nHalo:',nHalo)
    progress = -1
    
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs: ', progress, '%')
        nsub[i] = len([item for item in f[str(i)].values() if isinstance(item, h5py.Group)])

    print('nsub:',nsub)
    print('len(nsub):',len(nsub))

    bins = np.amax(nsub)

    H, bin_edge = np.histogram(nsub, bins)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bin_edge[:-1], H+1, color='r')
    plt.title('The Number of 0.1 Halos in a 0.2 Halo')

    ax.grid(True)
    ax.set_yscale('log')
    ax.set_xlabel('Number of Sub-Halos')
    ax.set_ylabel('Frequency + 1')

    plt.savefig('Plots/NumberOccupancy.png', dpi=fig.dpi)

    plt.close()

    return

def MaxMassCompPlot(halo_snapshot):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    #print('f.keys():',list(f.keys()))
    datasets = [item for item in f.values() if isinstance(item, h5py.Dataset)]
    groups = [item for item in f.values() if isinstance(item, h5py.Group)]
    print('len(groups):',len(groups))
    print('len(datasets):',len(datasets))
    base_list=len(f['/Halo_IDs/'])
    nHalo = len(groups)
    print('nHalo:',nHalo)
    main_mass = np.empty(nHalo,dtype=np.int32)
    nHaloSub=np.empty(nHalo,dtype=np.int32)
    progress = -1

    # First loop to count halos
    iSub=0
    nSub=0
    for i in range(0,nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress main_IDs: ', progress, '%')
        main_mass[i]=f['/'+str(i)].attrs['halo_nPart']
        nHaloSub[i] = len([item for item in f['/'+str(i)].values() if isinstance(item, h5py.Group)])
    nSub=np.sum(nHaloSub)

    # Array to hold subhalo masses
    sub_mass=np.empty(nSub,dtype=np.int32)
    split_sub_mass = np.empty(nHalo,dtype=np.int32)
    tot_sub_mass=np.empty(nHalo,dtype=np.int32)
    iSub=0
    
    progress = -1
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs: ', progress, '%')
        for j in range(0,nHaloSub[i]):
            sub_mass[iSub]=f['/'+str(i)+'/'+str(j)].attrs['subhalo_nPart']
            #print('len(sub_mass):',len(sub_mass))
            #print('sub_mass:',sub_mass)
            iSub +=1

    print('nHaloSub:',nHaloSub)
    print('len(nHaloSub):',len(nHaloSub))
    split_sub_mass = np.split(sub_mass, np.cumsum(nHaloSub)[:-1])
    print('len(split_sub_mass):',len(split_sub_mass))
    #print('split_sub_mass:',split_sub_mass)
    max_sub_mass = [max(x, default=0) for x in split_sub_mass]
    Main_mass = [x if x in main_mass else 0 for x in max_sub_mass]
    #print('sub_mass:',sub_mass)
    #print('len(nHaloSub):',len(main_mass))
    #print('main_mass:',main_mass)
    print('len(main_mass):',len(main_mass))
    print('main_mass:',main_mass)
    print('len(sub_mass):',len(sub_mass))
    print('sub_mass:',sub_mass)
    #print('split_sub_mass:',split_sub_mass)
    #print('max_sub_mass:',(max_sub_mass))
    print('len(max_sub_mass):',len((max_sub_mass)))

    mass_largest = [x for x in max_sub_mass if x >=10000]
    print('mass_largest:', mass_largest)

    fig = plt.figure()
    ax = fig.add_subplot(111)
 
    ax.plot(main_mass, max_sub_mass,'.', color='g')

    ax.grid(True)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_title('The Largest 0.1 Halo')
    ax.set_xlabel(r'$M_{0.2}$')
    ax.set_ylabel(r'$M_{0.1}$')

    plt.savefig('Plots/MaxMassComp.png', dpi=fig.dpi)

    plt.clf()
    
    return
