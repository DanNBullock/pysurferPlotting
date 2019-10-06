#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:47:42 2017

@author: plab
"""

import os as os

import numpy as np
from surfer import io
from surfer import Brain
import mayavi
import nibabel as nib


import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



# the histogram of the data




def singleNiftiHeatmap(niftiPath1, subject_id, subjects_dir):   
        
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
        if averageXcoord>xDim/2:
            side='lh'
        else:
            side='rh'
            
        dirContent=os.listdir(subjects_dir+subject_id)
        fsindices = [i for i, s in enumerate(dirContent) if 'freesurfer' in s]
        fsDir=NiftiPreDir+dirContent[fsindices[0]]    
            
        brain = Brain(subject_id='output', hemi=side, surf="inflated", background="white", subjects_dir=fsDir)
        
        #set this up for each system, I guess
        os.environ['FREESURFER_HOME'] ='/Applications/freesurfer/'
        os.environ['SUBJECTS_DIR'] =fsDir
        #print os.environ['SUBJECTS_DIR']
        

   
        steps= [float(i) for i in [0,10,.1]]
        
        #fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=[4])
        fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id='output', projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps)


        #thresh=max(fib1_surf_data)*.05
        
        fib1Max=max(fib1_surf_data)
#        fib1Dev=np.std(fib1_surf_data[(fib1_surf_data!=0)])
#        fib2Dev=np.std(fib2_surf_data[(fib2_surf_data!=0)])
#        fib1Ratio=fib1Max/fib1Dev
#        fib2Ratio=fib2Max/fib2Dev
#        fib1mean=np.true_divide(fib1_surf_data.sum(),(fib1_surf_data!=0).sum())
#        fib2mean=np.true_divide(fib2_surf_data.sum(),(fib2_surf_data!=0).sum())
        
        dataSort=fib1_surf_data
        dataSort=dataSort[dataSort>0]
        print dataSort
        
    
         
        #thresh=dataSort[twentiethPrctIndx]
        thresh=5
        print thresh


        

          

        
        brain.add_data(fib1_surf_data, 0, fib1Max ,.333,colormap="gist_rainbow", alpha=.6,
                       hemi=side)
                

        figdir=os.path.join(subjects_dir,subject_id ,'Figures/' )
        filename1=os.path.basename(niftiPath1)
        #i have no idea why this picks the first ., but it works.
        filedotIndex1=filename1.index('.')
        filetitle1=filename1[0:filedotIndex1]
        
        if not os.path.exists(figdir):
            os.makedirs(figdir)
            #views=[(-45,27.5),(-45,90),(45,90),(-90,135),(-225,27.5),(-315,27.5),(-45,-27.5),(-135,-27.5),(-225,-27.5),(-315,-27.5),(0,0),(0,-27.5),(0,27.5)]
        if side=='rh':
            views=[(-45,27.5),(-45,110),(45,90),(-90,135),(-315,27.5),(90,-180),(0,27.5),(-100,80),(-120,105),(0,90),(90,90)]
        else:
            views=[(-135,27.5),(-135,110),(135,90),(90,-135),(315,-27.5),(90,-180),(0,-27.5),(100,-80),(120,-105),(0,-90),(90,90)]
        
        for iviews in np.arange(len(views)):

            mayavi.mlab.view(azimuth=views[iviews][0], elevation=views[iviews][1])
            figName=filetitle1+'_heatmap'+str(iviews)+'.png'
            curFileName=os.path.join(figdir ,figName)
            brain.save_single_image(curFileName)
            
            
            
            




#niftis=[6, 7, 8, 9, 12, 13]

# 
subjIDs=['sub-103','sub-207','sub-304']
#subjIDs=['105115']
#NiftiDir='/Users/plab/Downloads/EnsOverlap2/105115/105115/'
#subject_id='105115'
#subjects_dir='/Users/plab/Downloads/EnsOverlap/'
#regPath='/Users/plab/Downloads/EnsOverlap/105115/mri/transforms/reg.mni152.1mm.dat'
#    NiftiDir='/Users/plab/Downloads/EnsOverlap2/DataandDerivatives/'+subjIDs[iG]+'/EndpointVolumes/'

for iG in np.arange(np.shape(subjIDs)[0]):
    subject_id=subjIDs[iG]
    subjects_dir='/Users/plab/Downloads/SophiaProject/proj-5a74d1e26ed91402ce400cca/'
    
    
    NiftiPreDir=subjects_dir+subject_id+'/'
    
    dirContent=os.listdir(NiftiPreDir)
    
    fsindices = [i for i, s in enumerate(dirContent) if 'freesurfer' in s]
    fsDir=NiftiPreDir+dirContent[fsindices[0]]+'/output/'
    
    regPath=fsDir+'mri/transforms/reg.mni152.1mm.dat'
    
    indices = [i for i, s in enumerate(dirContent) if 'rois' in s]
    
    NiftiDir=NiftiPreDir+dirContent[indices[0]]+'/'
    files = []
    for r, d, f in os.walk(NiftiDir):
        for file in f:
            if '.nii.gz' in file:
                files.append(os.path.join(r, file))

    niftis=np.arange(np.shape(files)[0])
    comparisonDim=np.shape(niftis)


    
    for ii in [0]:
        for ib in np.arange(comparisonDim[0]):
    
            singleNiftiHeatmap(files[ib], subject_id, subjects_dir)
    