import numpy as np
import matplotlib.pyplot as plt
import h5py

def SubHaloPosVel(halo_snapshot, occupancy_l, occupancy_h):
    f=h5py.File(halo_snapshot, 'r', driver='core')
    nHalo = len([item for item in f.values() if isinstance(item, h5py.Group)])
    nHaloSub=np.empty(nHalo,dtype=np.int32)
    Halo_Pos =  np.empty(nHalo,dtype=np.float32)
    halo_pos =  np.empty(nHalo,dtype=np.float32)
    Halo_Vel =  np.empty(nHalo,dtype=np.float32)
    halo_vel=  np.empty(nHalo,dtype=np.float32)
    Sigma =  np.empty(nHalo,dtype=np.float32)

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
        Vel_Halo=f['/'+str(i)+'/'+'Halo_Vel'].value
        Vel_Halo_Mean=f['/'+str(i)+'/'+'mean_vel'].value
        halo_vel =  np.empty(len(Vel_Halo),dtype=np.float32)
        for j in range (0, len(Pos_Halo)):
            halo_pos[j]= np.sqrt(np.power(Pos_Halo[j,0],2,dtype=np.float32) + np.power(Pos_Halo[j,1],2,dtype=np.float32) +np.power(Pos_Halo[j,2],2,dtype=np.float32), dtype=np.float32)
        max_halo_pos = np.max(halo_pos)
        Halo_Pos[i]= max_halo_pos
        for j in range (0, len(Vel_Halo)):
            halo_vel[j]= np.sqrt(np.power(Vel_Halo[j,0],2,dtype=np.float32) + np.power(Vel_Halo[j,1],2,dtype=np.float32) +np.power(Vel_Halo[j,2],2,dtype=np.float32), dtype=np.float32)
        mean_halo_vel = np.mean(halo_vel)
        Halo_Vel[i]= mean_halo_vel
        U = Vel_Halo - Vel_Halo_Mean
        Sigma[i] = np.dot(U[j,:],U[j,:])
        if occupancy_l <= nhalosub and nhalosub <= occupancy_h:
            nHaloSub[i] = nhalosub 
    nSub=np.sum(nHaloSub)

    Sub_Halo_Pos = np.empty(nSub,dtype=np.float32)
    Sub_Halo_Vel = np.empty(nSub,dtype=np.float32)
    RatioPos = np.empty(nSub,dtype=np.float32)
    RatioVel = np.empty(nSub,dtype=np.float32)
        
    progress = -1
    for i in range(0, nHalo):
        previous_progress = progress
        progress = int(i/nHalo * 100)
        if progress != previous_progress:
            print('Progress sub_IDs: ', progress, '%')
        for j in range(0,nHaloSub[i]):
            Pos_Sub=f['/'+str(i)+'/'+str(j)+'/subhalo_mean_pos'].value
            Sub_Halo_Pos = np.sqrt(np.power(Pos_Sub[0],2,dtype=np.float32) + np.power(Pos_Sub[1],2,dtype=np.float32) +np.power(Pos_Sub[2],2,dtype=np.float32), dtype=np.float32)
            Vel_Sub=f['/'+str(i)+'/'+str(j)+'/subhalo_mean_vel'].value
            Sub_Halo_Vel = np.sqrt(np.power(Vel_Sub[0],2,dtype=np.float32) + np.power(Vel_Sub[1],2,dtype=np.float32) +np.power(Vel_Sub[2],2,dtype=np.float32), dtype=np.float32)
            RatioPos[iSub] = Sub_Halo_Pos/Halo_Pos[i]
            RatioVel[iSub] = abs(Sub_Halo_Vel-Halo_Vel[i])/np.sqrt(Sigma[i])
            iSub += 1
       
    RatioPos =  RatioPos[~np.isnan(RatioPos)]
    RatioVel =  RatioVel[~np.isnan(RatioVel)]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(RatioPos, RatioVel, color='r', marker='.')
    ax.grid(True)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(10**-4, 0.2*10**2)
    
    ax.set_xlabel('Sub-halo position in Main halo')
    ax.set_ylabel('Sub-halo relative speed')

    plt.savefig('Plots/SubHaloPosVel_Occ='+str(occupancy_l)+'-'+str(occupancy_h)+'.png', dpi=fig.dpi)

    plt.clf()
    
    return
