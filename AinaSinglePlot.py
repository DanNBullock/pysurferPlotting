#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:34:20 2017

@author: plab
"""

import numpy as np
from surfer import io
from surfer import Brain
import mayavi
import nibabel as nib
import os
import pandas as pd


#NiftiDir='/Users/plab/Downloads/EnsOverlap2/105115/105115/'
#subject_id='105115'
#subjects_dir='/Users/plab/Downloads/EnsOverlap/'
#regPath='/Users/plab/Downloads/EnsOverlap/105115/mri/transforms/reg.mni152.1mm.dat'

def readcsv(filename):
    data = pd.read_csv(filename,header=None) #Please add four spaces here before this line
    return(np.array(data))

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    minIndex=np.argmin(dist_2)
    return nodes[minIndex,:]

def nodes_in_range(node, nodes,maxRange):
    nodes = np.asarray(nodes)
    dist_2 = np.sqrt(np.sum((nodes - node)**2, axis=1))
    maxDistBool=dist_2<maxRange
    nodeIndexes=np.where(maxDistBool)
    return nodes[nodeIndexes[0],:]

def interpolate_coords_to_GM(fsPath,coords):
    ribbonVol=nib.load(fsPath+'mri/ribbon.mgz')
    ribbonData=ribbonVol.get_data()
    nifti_affine=ribbonVol.affine
    
    visualNifti=nib.load(os.environ['FREESURFER_HOME']+'average/mni305.cor.mgz')
    visual_affine=visualNifti.affine
    
    #set kernel here
    maxDist=8
    # grey version
    #leftBool=ribbonData==3
    #rightBool=ribbonData==42
    #white version
    leftBool=ribbonData==2
    rightBool=ribbonData==41
    leftIndx=np.where(leftBool)
    rightIndx=np.where(rightBool)
    greyIndx=np.concatenate((leftIndx,rightIndx),axis=1)
    non_zero_indexes_array=np.asarray(greyIndx)
    #does the array need to be transposed?
    non_zero_indexes_array_transposed=np.transpose(non_zero_indexes_array)
    affinedCoords=nib.affines.apply_affine(nifti_affine,non_zero_indexes_array_transposed)
    
    #cord size
    coordSize=coords.shape
    interpCoord=np.zeros(coordSize)
    for iInputCoords in range(coordSize[0]):
        candidateNodes=nodes_in_range(coords[iInputCoords], affinedCoords,maxDist)
        distArray=affinedCoords-coords[0]
        distBool=np.zeros(( len(distArray),1), dtype=bool)
        for iCoords in range(len(distArray)):
            distBool[iCoords]=np.sqrt(distArray[iCoords,0]**2+distArray[iCoords,1]**2+distArray[iCoords,2]**2)<maxDist
        nearOut=np.where(distBool)
        nearIndexes=nearOut[0]

        centroidVal= [np.mean(affinedCoords[nearIndexes,0]),np.mean(affinedCoords[nearIndexes,1]),np.mean(affinedCoords[nearIndexes,2])   ]
        centroidVal= np.mean(candidateNodes,axis=0)
        interpCoord[iInputCoords,:]=closest_node(centroidVal, candidateNodes)
        
    invert_nifti_affine=np.linalg.inv(nifti_affine)
    imgCoords=nib.affines.apply_affine(invert_nifti_affine,interpCoord)
    reconvertedCoords=nib.affines.apply_affine(visual_affine,imgCoords)
        
    return reconvertedCoords
    
def rewarp_coords(fsPath,coords):
    ribbonVol=nib.load(fsPath+'mri/ribbon.mgz')

    nifti_affine=ribbonVol.affine
    
    visualNifti=nib.load(os.environ['FREESURFER_HOME']+'average/mni305.cor.mgz')
    visual_affine=visualNifti.affine
    
   
        
    invert_nifti_affine=np.linalg.inv(nifti_affine)
    imgCoords=nib.affines.apply_affine(invert_nifti_affine,coords)
    reconvertedCoords=nib.affines.apply_affine(visual_affine,imgCoords)
        
    return reconvertedCoords
             
    


def MultiNiftiOverlap(niftiPath1, coordsPath ,subject_id, subjects_dir):   
        
        #for whatever reason, it seems that left and right are switched here.
        # at least this is what I infer from the nifti.  Right side niftis
        # have lower x values than left side niftis, which one wouldnt usualy
        # expect, but whatever.  It serves as the appropriate check
        niftiVol=nib.load(niftiPath1)
        niftiData=niftiVol.get_data()
        nonZeros=np.nonzero(niftiData)
        niftiDimensions=niftiData.shape
        xDim=niftiDimensions[1]
        averageXcoord=np.mean(nonZeros[0])
        if averageXcoord<xDim/2:
            side='lh'
        else:
            side='rh'
            
        print(niftiPath1)

            
        brain = Brain(subject_id=subject_id, hemi=side, surf="inflated", background="white", subjects_dir=subjects_dir)
        
        #set this up for each system, I guess
        os.environ['FREESURFER_HOME'] ='/Applications/freesurfer/'
        os.environ['SUBJECTS_DIR'] =subjects_dir
        #print os.environ['SUBJECTS_DIR']
        
        steps= [float(i) for i in [0,8,.1]]
        
        #NOTE:  FIBER 2 IS ACTUALLY THE ECOG PLOT HERE.
        
        
        
        
        epi_img = nib.load(niftiPath1)
        epi_img_data = epi_img.get_fdata()
        np.unique(epi_img_data.data)
        
        #unique, counts = np.unique(epi_img_data.data, return_counts=True)
        #dict(zip(unique, counts))
        
  
        fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps)
        #fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5)
     
        #Just scale the data I guess?  Something is being messed with in the data.

        #fib1Max=max(fib1_surf_data)


        #thresh1=np.median(fib1_surf_data[np.nonzero(fib1_surf_data)])
        #print(thresh1)
        #coords2=[]
        coords=readcsv(coordsPath)
        #padCoord=np.zeros((1,3))
        #coords2=np.concatenate((coords,padCoord),axis=0)
        
        
        fsPath=subjects_dir+subject_id+'/'
        interpolatedCoords=interpolate_coords_to_GM(fsPath,coords)
        #interpolatedCoords=rewarp_coords(fsPath,coords)
        
        

        brain.add_data(fib1_surf_data, 0, 1000  , 150 ,colormap="hsv", alpha=.68,
                       hemi=side)
        
        
#        colorwords=['red',
#                     'blue',
#                     'green',
#                     'black',
#                     'purple',
#                     'cyan',
#                     'yellow',
#                     'pink'
#                     'orange']
        
        
        
#        
#        coordsSize=interpolatedCoords.shape
#        
#        print(interpolatedCoords)
        
        
#        for iCoords in range(coordsSize[0]):
        brain.add_foci(interpolatedCoords, map_surface='white', color="black", scale_factor=.5)
#            print(interpolatedCoords[iCoords])
#            print(colorwords[iCoords])
#            brain.add_foci(interpolatedCoords[iCoords], map_surface='white', color=colorwords[iCoords], scale_factor=.5)


        figdir=os.path.join(subjects_dir,subject_id ,'Figures/' )
        filename1=os.path.basename(niftiPath1)

        #i have no idea why this picks the first ., but it works.
        filedotIndex1=filename1.index('.')
        
        coordName=os.path.basename(coordsPath)
        coordName=coordName.replace('.csv','')

        filetitle1=filename1[0:filedotIndex1]+'_'+coordName

        
        if not os.path.exists(figdir):
            os.makedirs(figdir)
            #views=[(-45,27.5),(-45,90),(45,90),(-90,135),(-225,27.5),(-315,27.5),(-45,-27.5),(-135,-27.5),(-225,-27.5),(-315,-27.5),(0,0),(0,-27.5),(0,27.5)]
        if side=='rh':
            views=[(-45,27.5),(-45,110),(45,90),(-90,135),(-315,27.5),(90,-180),(0,27.5),(-100,80),(-120,105),(0,90),(90,90)]
        else:
            views=[(-135,27.5),(-135,110),(135,90),(90,-135),(315,-27.5),(90,-180),(0,-27.5),(100,-80),(120,-105),(0,-90),(90,90)]
        
        for iviews in np.arange(len(views)):
            

            mayavi.mlab.view(azimuth=views[iviews][0], elevation=views[iviews][1])
            figName=filetitle1+ str(iviews)+'.png'
            curFileName=os.path.join(figdir ,figName)
            brain.save_single_image(curFileName)





ainaCoords=['/Users/plab/Downloads/Coords/R_FG_A.csv',
'/Users/plab/Downloads/Coords/R_FG_B.csv',
'/Users/plab/Downloads/Coords/R_IOS_A.csv',
'/Users/plab/Downloads/Coords/R_IOS_b.csv',
'/Users/plab/Downloads/Coords/R_ITS_A.csv',
'/Users/plab/Downloads/Coords/R_ITS_B.csv',
'/Users/plab/Downloads/Coords/R_STS_A.csv',
'/Users/plab/Downloads/Coords/R_STS_B.csv']

niftiDir='/Users/plab/Downloads/HCPAmalgumNiftis'

niftiDirContent=os.listdir(niftiDir)
niftiDirContent.remove('.DS_Store')


niftiDirData=pd.Series(niftiDirContent)
niftiBool=niftiDirData.str.contains('right', regex=False)


#subjIDs=['105115','110411','111312','113619','FP','HT','KK','MP']

#NiftiDir='/Users/plab/Downloads/EnsOverlap2/105115/105115/'
#subjIDs='105115'
#subjects_dir='/Users/plab/Downloads/EnsOverlap/'
#regPath='/Users/plab/Downloads/EnsOverlap/105115/mri/transforms/reg.mni152.1mm.dat'
sourceNifti='/Users/plab/Downloads/Ecog Niftis-selected/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c.nii'

#NiftiDir='/Users/plab/Downloads/EnsOverlap2/DataandDerivatives/ICBM2009c_asym_nlin/EndpointVolumes/'
#AinaDir=''
subject_id='ICBM2009c_asym_nlin'
subjects_dir='/Users/plab/Downloads/EnsOverlap/'
regPath='/Users/plab/Downloads/EnsOverlap/ICBM2009c_asym_nlin/mri/transforms/reg.mni152.1mm.dat'
    
#coord_csv_to_sphere_roi_nifti('/Users/plab/Desktop/PythonCode/Coords/R_FG_A.txt',5,sourceNifti)

#for iTracts in range(len(tractList)/2,len(tractList)):
for iEcog in range(len(ainaCoords)):

     curCoords=ainaCoords[iEcog]
     
     for iNifti in range(len(niftiDirContent)):
         
         curNifti=niftiDir+'/'+niftiDirContent[iNifti]
         MultiNiftiOverlap(curNifti, curCoords, subject_id, subjects_dir)
    
    
  
    
    
    