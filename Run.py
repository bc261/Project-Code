import mainBranchLength as mBL
import forestgrowth as fg
import ProgDescHistogram as PDH
import logMgrowth as lMg
import massfluctuation as mf
import MassComp as MC
import SubHaloPos as SHP
import SubHaloVel as SHV
import SubPosVel as SPV
import pickle
import numpy

#Linking length = 0.1
#with open('treedata.pck', 'rb') as pfile1:
#    tree_data_1 = pickle.load(pfile1)

#Linking length = 0.2
#with open('treedata10.pck', 'rb') as pfile2:
#    tree_data_2 = pickle.load(pfile2)

#with open('snapdatatime.pck', 'rb') as pfile3:
#    snap_data = pickle.load(pfile3)


#lMg.mainLogM(tree_data, snap_data, cutoff=10)

#mf.mainmfluc(tree_data, snap_data, cutoff=1000)

#mBL.mainBranchLengthCompPlot(tree_data, cutoff=None)

#PDH.directprogdeschist(tree_data, cutoff=1000)

#PDH.progdescsum(tree_data_1, tree_data_2, cutoff=1000)

#PDH.progdescstat(tree_data_1, tree_data_2, cutoff=1000)

#PDH.progdeschistComp(tree_data_1, tree_data_2, cutoff=1000)

#lMg.mainLogMComp(tree_data_1, tree_data_2, snap_data, cutoff=1000)

#mf.mainmflucComp (tree_data_1, tree_data_2, snap_data, cutoff=1000)

#MC.SumMassCompPlot(halo_snapshot='halos_061.hdf5')

#SHP.SubHaloPos(halo_snapshot= 'halos_061.hdf5', occupancy_l=2, occupancy_h=4)

#SHP.SubHaloPosComp(halo_snapshot='halos_061.hdf5', occupancy1_l=0, occupancy1_h=1, occupancy2_l=2, occupancy2_h=4)

#SHV.SubHaloVel(halo_snapshot= 'halos_061.hdf5', occupancy_l=0, occupancy_h=1)

#SHV.SubHaloVelComp(halo_snapshot= 'halos_061.hdf5', occupancy1_l=0, occupancy1_h=1, occupancy2_l=2, occupancy2_h=4)

#MC.NumOcup(halo_snapshot='halos_061.hdf5')

#MC.MaxMassCompPlot(halo_snapshot='halos_061.hdf5')

SPV.SubHaloPosVel(halo_snapshot='halos_061.hdf5', occupancy_l=0, occupancy_h=1)
