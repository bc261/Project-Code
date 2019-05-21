import numpy as np
import matplotlib.pyplot as plt
import h5py

def SubHaloPos(halo_snapshot, occupancy_l, occupancy_h):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    nHalo = len([item for item in f.values() if isinstance(item, h5py.Group)])
    nHaloSub=np.empty(nHalo,dtype=np.int32)
    Halo_Pos =  np.empty(nHalo,dtype=np.float32)
    halo_pos =  np.empty(nHalo,dtype=np.float32)
    
    progress = -1
    iSub=0
    nSub=0
    for i in range(0,nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress main_IDs: ', progress, '%')
        nhalosub = len([item for item in f['/'+str(i)].values() if isinstance(item, h5py.Group)])
        Pos_Halo=f['/'+str(i)+'/'+'Halo_Pos'].value
        halo_pos =  np.empty(len(Pos_Halo),dtype=np.float32)
        for j in range (0, len(Pos_Halo)):
            halo_pos[j]= np.sqrt(np.power(Pos_Halo[j,0],2,dtype=np.float32) + np.power(Pos_Halo[j,1],2,dtype=np.float32) +np.power(Pos_Halo[j,2],2,dtype=np.float32), dtype=np.float32)
        max_halo_pos = np.max(halo_pos)
        Halo_Pos[i]= max_halo_pos
        if occupancy_l <= nhalosub and nhalosub <= occupancy_h:
            nHaloSub[i] = nhalosub 
    nSub=np.sum(nHaloSub)

    Sub_Halo_Pos = np.empty(nSub,dtype=np.float32)
    Ratio = np.empty(nSub,dtype=np.float32)
     
    progress = -1
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs: ', progress, '%')
        for j in range(0,nHaloSub[i]):
            Pos_Sub=f['/'+str(i)+'/'+str(j)+'/subhalo_mean_pos'].value
            Sub_Halo_Pos = np.sqrt(np.power(Pos_Sub[0],2,dtype=np.float32) + np.power(Pos_Sub[1],2,dtype=np.float32) +np.power(Pos_Sub[2],2,dtype=np.float32), dtype=np.float32)
            Ratio[iSub] = Sub_Halo_Pos/Halo_Pos[i]
            iSub += 1
  
    Ratio =  Ratio[~np.isnan(Ratio)]
    
    bins = 50

    H , bin_edge = np.histogram(Ratio, bins, density=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bin_edge[:-1], H , color='r')
    ax.grid(True)
    ax.set_xlim(0,1)
   
    plt.title('Positions of Sub-halos of occupancy '+ str(occupancy_l)+'-'+ str(occupancy_h)+ ' compared to the Main Halos centre')
    ax.set_xlabel('Sub-Halo position from Main Halo centre divide by radius of Main Halo ')
    ax.set_ylabel('Density')

    plt.savefig('Plots/SubHaloPos_Occ='+str(occupancy_l)+'-'+str(occupancy_h)+'.png', dpi=fig.dpi)

    plt.clf()
    
    return

def SubHaloPosComp (halo_snapshot, occupancy1_l, occupancy1_h, occupancy2_l, occupancy2_h):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    nHalo = len([item for item in f.values() if isinstance(item, h5py.Group)])
    nHaloSub_1=np.empty(nHalo,dtype=np.int32)
    nHaloSub_2=np.empty(nHalo,dtype=np.int32)
    Halo_Pos =  np.empty(nHalo,dtype=np.float32)
    halo_pos =  np.empty(nHalo,dtype=np.float32)

    progress = -1
    nSub_1=0
    nSub_2=0
    for i in range(0,nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress main_IDs: ', progress, '%')
        nhalosub = len([item for item in f['/'+str(i)].values() if isinstance(item, h5py.Group)])
        Pos_Halo=f['/'+str(i)+'/'+'Halo_Pos'].value
        halo_pos =  np.empty(len(Pos_Halo),dtype=np.float32)
        for j in range (0, len(Pos_Halo)):
            halo_pos[j]= np.sqrt(np.power(Pos_Halo[j,0],2,dtype=np.float32) + np.power(Pos_Halo[j,1],2,dtype=np.float32) +np.power(Pos_Halo[j,2],2,dtype=np.float32), dtype=np.float32)
        max_halo_pos = np.max(halo_pos)
        Halo_Pos[i]= max_halo_pos
        if occupancy1_l <= nhalosub and nhalosub <= occupancy1_h:
            nHaloSub_1[i] = nhalosub
        if occupancy2_l <= nhalosub and nhalosub <= occupancy2_h:
            nHaloSub_2[i] = nhalosub
            
    nSub_1=np.sum(nHaloSub_1)
    nSub_2=np.sum(nHaloSub_2)

    Sub_Halo_Pos_1 = np.empty(nSub_1,dtype=np.float32)
    Sub_Halo_Pos_2 = np.empty(nSub_2,dtype=np.float32)
    Ratio_1 = np.empty(nSub_1,dtype=np.float32)
    Ratio_2 = np.empty(nSub_2,dtype=np.float32)

     
    progress = -1
    iSub=0
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs for first halo: ', progress, '%')
        for j in range(0,nHaloSub_1[i]):
            Pos_Sub_1=f['/'+str(i)+'/'+str(j)+'/subhalo_mean_pos'].value
            Sub_Halo_Pos_1 = np.sqrt(np.power(Pos_Sub_1[0],2,dtype=np.float32) + np.power(Pos_Sub_1[1],2,dtype=np.float32) +np.power(Pos_Sub_1[2],2,dtype=np.float32), dtype=np.float32)
            Ratio_1[iSub] = Sub_Halo_Pos_1/Halo_Pos[i]
            iSub += 1

    progress = -1
    iSub=0
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs for second halo: ', progress, '%')
        for j in range(0,nHaloSub_2[i]):
            Pos_Sub_2=f['/'+str(i)+'/'+str(j)+'/subhalo_mean_pos'].value
            Sub_Halo_Pos_2 = np.sqrt(np.power(Pos_Sub_2[0],2,dtype=np.float32) + np.power(Pos_Sub_2[1],2,dtype=np.float32) +np.power(Pos_Sub_2[2],2,dtype=np.float32), dtype=np.float32)
            Ratio_2[iSub] = Sub_Halo_Pos_2/Halo_Pos[i]
            iSub += 1
        
  
    Ratio_1 =  Ratio_1[~np.isnan(Ratio_1)]
    Ratio_2 =  Ratio_2[~np.isnan(Ratio_2)]

    bins = 50

    H_1 , bin_edge_1 = np.histogram(Ratio_1, bins, density=True)
    H_2 , bin_edge_2 = np.histogram(Ratio_2, bins, density=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bin_edge_1[:-1], H_1 , color='r', label='occupancy '+str(occupancy1_l)+'-'+str(occupancy1_h))
    ax.plot(bin_edge_2[:-1], H_2 , color='g', label='occupancy '+str(occupancy2_l)+'-'+str(occupancy2_h))
    ax.grid(True)
    ax.set_xlim(0,1)
   
    plt.title('Positions of Sub-halos of different occupancy compared to the Main Halos centre')
    ax.set_xlabel('Sub-Halo position from Main Halo centre divide by radius of Main Halo ')
    ax.set_ylabel('Density')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    plt.savefig('Plots/SubHaloPosComp.png', dpi=fig.dpi)

    plt.clf()
    
    return
