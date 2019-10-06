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

def MultiNiftiOverlap(niftiPath1, niftiPath2, subject_id, subjects_dir):   
        
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
            
        brain = Brain(subject_id=subject_id, hemi=side, surf="white", background="white", subjects_dir=subjects_dir)
        
        #set this up for each system, I guess
        os.environ['FREESURFER_HOME'] ='/Applications/freesurfer/'
        os.environ['SUBJECTS_DIR'] =subjects_dir
        #print os.environ['SUBJECTS_DIR']
        
        
        #NOTE: There is also another optional variable that you can put in here, projarg=steps where steps= [float(i) for i in [0,8,.1]]
        # what this will do is specify how deep into the cortex it is going to look for data.  If your plots are coming out kind of sparse, feel free to try this
        fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5)
        fib2_surf_data = io.project_volume_data(filepath=niftiPath2, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5)
     
        #Here we are setting the maximum color value for the color map
        fib1Max=max(fib1_surf_data)
        fib2Max=max(fib2_surf_data)
        
        #here we are setting the lower bound threshold for this plotting instance.  Anything lower than this is zeroed out on the plot
        thresh1=fib1Max*.005
        thresh2=fib2Max*.005
        
        #note the value that preceeds the fiber max.  This stretches the color map so that more striking colors are used
        brain.add_data(fib1_surf_data, -fib1Max*.75, fib1Max  ,thresh2 ,colormap="Reds", alpha=.5,
                       hemi=side)         
        
        brain.add_data(fib2_surf_data, -fib2Max*.75 , fib2Max , thresh1, colormap="Blues", alpha=.5,
                       hemi=side)


        #sets the figure directory
        figdir=os.path.join(subjects_dir,subject_id ,'Figures/' )
        filename1=os.path.basename(niftiPath1)
        filename2=os.path.basename(niftiPath2)
        #i have no idea why this picks the first ., but it works.
        filedotIndex1=filename1.index('.')
        filedotIndex2=filename2.index('.')
        filetitle1=filename1[0:filedotIndex1]
        filetitle2=filename2[0:filedotIndex2]
        
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        
        #changes the views depending on the hemisphere used
        if side=='rh':
            views=[(-45,27.5),(-45,110),(45,90),(-90,135),(-315,27.5),(90,-180),(0,27.5),(-100,80),(-120,105),(0,90),(90,90)]
        else:
            views=[(-135,27.5),(-135,110),(135,90),(90,-135),(315,-27.5),(90,-180),(0,-27.5),(100,-80),(120,-105),(0,-90),(90,90)]
    
        #iterates through views and saves down
        for iviews in np.arange(len(views)):
            #this mayavi bit is why this won't work on clusters.  It needs an interactive window
            #to manipulate in order to generate and iterate through images
            mayavi.mlab.view(azimuth=views[iviews][0], elevation=views[iviews][1])
            figName=filetitle1+'_to_'+filetitle2+ str(iviews)+'.png'
            curFileName=os.path.join(figdir ,figName)
            brain.save_single_image(curFileName)


# this is an inelegant way of making a list of the relevant files you are going to iterate through.
#  A better way might be to use os.listdir to get a list of all of the contents 
#  Replace these with the names of the niftis that are relevant to your project
nameList=['amalgum_Left_Arcuate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Arcuate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Baum_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Baum_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Cingulum_Cingulate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Cingulum_Cingulate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Cingulum_Hippocampus_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Cingulum_Hippocampus_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Corticospinal_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Corticospinal_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_IFOF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_IFOF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_ILF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_ILF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_MdLF-SPL_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_MdLF-SPL_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_MdLF-Ang_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_MdLF-Ang_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Meyer_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Meyer_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_pArc_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_pArc_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_SLF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_SLF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Thalamic_Radiation_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Thalamic_Radiation_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_TPC_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_TPC_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_Uncinate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_Uncinate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Left_VOF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Left_VOF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Arcuate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Arcuate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Baum_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Baum_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Cingulum_Cingulate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Cingulum_Cingulate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Cingulum_Hippocampus_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Cingulum_Hippocampus_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Corticospinal_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Corticospinal_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_IFOF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_IFOF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_ILF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_ILF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_MdLF-SPL_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_MdLF-SPL_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_MdLF-Ang_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_MdLF-Ang_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Meyer_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Meyer_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_pArc_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_pArc_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_SLF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_SLF_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Thalamic_Radiation_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Thalamic_Radiation_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_TPC_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_TPC_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_Uncinate_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_Uncinate_gaussian_15mm_RAS_FiberEndpoint.nii.gz',
'amalgum_Right_VOF_gaussian_15mm_LPI_FiberEndpoint.nii.gz',
'amalgum_Right_VOF_gaussian_15mm_RAS_FiberEndpoint.nii.gz']

# here you list the indexes of the file names in the nameList that you would like to see plotted together.  Be sure to remember 0 indexing and Left/Right issues.
comparisons=[[47,59],[15,27],[46,58],[14,26], [19,25], [20,26], [51,57],[52,58],[15,13],[16,14],[47,45],[48,46],[13,25],[14,26],[45,57],[46,58],[15,19],[16,20],[47,51],[48,52] ]
comparisonDim=np.shape(comparisons) 

#if you want to iterate through multiple subjects you can list them here.
#unfortuantely for us pysurfer makes some assumptions about the directory structure of 
#our freesurfer output, and presumes that there is an omnibus freesurfer directory.
#As a consequence of this you may need to change the subjectID to "output" and go one subject at a time.
# to iterate through subjects you would then alter the subjects_dir input variable
subjIDs=['105115']
#some example paths
#NiftiDir='/Users/plab/Downloads/EnsOverlap2/105115/105115/'
#subjIDs='105115'
#subjects_dir='/Users/plab/Downloads/EnsOverlap/'
#regPath='/Users/plab/Downloads/EnsOverlap/105115/mri/transforms/reg.mni152.1mm.dat'

for iG in np.arange(np.shape(subjIDs)[0]):
    #path to this subject's endpoint density volumes
    NiftiDir='/Users/plab/Downloads/EnsOverlap2/DataandDerivatives/'+subjIDs[iG]+'/EndpointVolumes/'
    #freesurfer subject ID.  Unfortunately, due to how brainlife structures output, all of our subject's names are "output"
    subject_id=subjIDs[iG]
    #Path to freesurfer subject directory.  Will output figures here and set environmental variables accordingly
    subjects_dir='/Users/plab/Downloads/EnsOverlap/'
    #path to registration file
    regPath='/Users/plab/Downloads/EnsOverlap/'+subjIDs[iG]+'/mri/transforms/reg.mni152.1mm.dat'
    
    for ii in [0]:
        for ib in np.arange(comparisonDim[0]):
    
            MultiNiftiOverlap(NiftiDir+nameList[comparisons[ib][0]+ii], NiftiDir+nameList[comparisons[ib][1]+ii], subject_id, subjects_dir)
    
    
  
    
    
    