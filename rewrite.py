#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:01:51 2022

@author: dan
"""

import numpy as np
from surfer import io
from surfer import Brain
import mayavi
import nibabel as nib

testNifti='/home/dan/Desktop/accSeg2/segViaDocker/output/seg/subCompTiles/Brain-Stem_TO_ctx-lh-superiorfrontal.nii.gz'
niftiPath1=testNifti
testFSDir='/home/dan/Desktop/accSeg2/testdata/output'
targetSubj_FS_Dir=testFSDir

def singleNiftiHeatmap_clean(niftiPath1,targetSubj_FS_Dir,surfType='inflated', outFigDir=None,figName=None):
    import nibabel as nib
    from nilearn.image import crop_img, resample_img
    import os
    import dipy.segment.mask as mask
    fsPathOut=os.path.split(targetSubj_FS_Dir)
    #not really, but for our purposes yes
    subjects_dir=fsPathOut[0]
    #curr subject proxy
    subject_id=fsPathOut[1]
    
    if isinstance(niftiPath1,str):
        inputNifti=nib.load(niftiPath1)
    else:
        inputNifti=niftiPath1
    
    #load up a freesurferf image to drop cerebellar and spinal rois, which are not surafce plottable
    fsAtlas=nib.load(os.path.join(targetSubj_FS_Dir,'mri/aparc.DKTatlas+aseg.nii.gz'))
    
    relabeledAtlas=wmaPyTools.roiTools.changeMultpleLabelValues(fsAtlas,[7,46,8,47,16],123456)
    
    #mask out the input nifti to find where the relevant data boundaries are
    reshapedAtlas=resample_img(relabeledAtlas,target_affine=inputNifti.affine,target_shape=inputNifti.shape)
      
    #get input image data
    inputNiftiData=inputNifti.get_data()
    
    #mask out the out of mask data
    inputNiftiData[reshapedAtlas==123456]=0
    
    #select the voxels from the input nifti that meet this criterion
    fsMaskedData=nib.Nifti1Image(inputNiftiData, inputNifti.affine,inputNifti.header)
    
    #first, remove any voxels from the input that are not plotted to surface,
    #at least for the purposes of hemisphere inference.
    subjectSpaceNiftiBounds=wmaPyTools.roiTools.subjectSpaceMaskBoundaryCoords(fsMaskedData)
    
    if np.all(subjectSpaceNiftiBounds[:,0]>0):
        side='rh'
        print('Right hemisphere nifti map detected')
    elif np.all(subjectSpaceNiftiBounds[:,0]<0):
        side='lh'
        print('Left hemisphere nifti map detected')
    else:
        side='split'
        print('Left AND right hemisphere nifti map detected, "both" will be plotted')
        print('If this was not the intended case, clean input nifti and resubmit')
        
    
    brain = Brain(subject_id=subject_id, hemi=side, surf=surfType, background="white", subjects_dir=subjects_dir)
   
    #not necessary?
    #set this up for each system, I guess
    #os.environ['FREESURFER_HOME'] ='/Applications/freesurfer/'
    #a workaround that may or may not work
    os.environ['SUBJECTS_DIR'] =subjects_dir
  
    steps= [float(i) for i in [0,10,.1]]
    
    thresh=.333
    
    #fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=[4])
    #yes, this is a mess
    if side in ['both','split']:
        for iSides in ['lh','rh']:
            #tempHemi=iSides
            if isinstance(niftiPath1, str):
                #for whatever reason it has to be a file on disk
                fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=iSides,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps,verbose=True)
            else:
                print ('freesurfer stubbornly requires input nifti to be on disk, using \n ' +inputNifti.get_filename() + '\n as path')
                fib1_surf_data = io.project_volume_data(filepath=inputNifti.get_filename(), hemi=iSides,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps,verbose=True)
            fib1Max=max(fib1_surf_data)

            brain.add_data(fib1_surf_data, 0, fib1Max ,thresh,colormap="gist_rainbow", alpha=.6, hemi=iSides)
    else:
        if isinstance(niftiPath1, str):
            #for whatever reason it has to be a file on disk
            fib1_surf_data = io.project_volume_data(filepath=niftiPath1, hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps,verbose=True)
        else:
            print ('freesurfer stubbornly requires input nifti to be on disk, using \n ' +inputNifti.get_filename() + '\n as path')
            fib1_surf_data = io.project_volume_data(filepath=inputNifti.get_filename(), hemi=side,  subject_id=subject_id, projsum='max', smooth_fwhm=5, projmeth='dist',  projarg=steps,verbose=True)
        fib1Max=max(fib1_surf_data)

        brain.add_data(fib1_surf_data, 0, fib1Max ,thresh,colormap="gist_rainbow", alpha=.6, hemi=side)
        

    
    if not outFigDir==None:
        if not os.path.exists(outFigDir):
            os.makedirs(outFigDir)

    #it's probably possible to use a clsutering method to establish clusters,
    #find their centroids, and then compute azimuth and elevation for these.
    if side=='rh':
        views=[(-45,27.5),(-45,110),(45,90),(-90,135),(-315,27.5),(90,-180),(0,27.5),(-100,80),(-120,105),(0,90),(90,90)]
    else:
        views=[(-135,27.5),(-135,110),(135,90),(90,-135),(315,-27.5),(90,-180),(0,-27.5),(100,-80),(120,-105),(0,-90),(90,90)]
    
    
    for iterator,iviews in enumerate(views):
        #this was the old way
        #mayavi.mlab.view(azimuth=views[iviews][0], elevation=views[iviews][1])
        view={'azimuth':iviews[0], 'elevation':iviews[1]}
        brain.show_view(view=view)
        #parse or create a fig name
        if figName==None:
            pathParts=os.path.split(inputNifti.get_filename())
            curFigStem=pathParts[1].replace('.nii.gz','')
        else:
            curFigStem=figName
            
        figName=curFigStem+'_heatmap'+str(iviews)+'.png'
        curFileName=os.path.join(outFigDir ,figName)
        brain.save_single_image(curFileName)
