import SimpleITK as sitk
import numpy as np
import time,shutil
from glob import iglob
import os
import argparse

def GetListFiles(folder_path, file_extension):
    """Return a list of files in folder_path finishing by file_extension"""
    file_list = []
    for extension_type in file_extension:
        file_list += search(folder_path,file_extension)[extension_type]
    return file_list

def GetPatients(folder_path, time_point='T1'):
    """Return a dictionary with patient id as key"""

    file_extension = ['.nii.gz','.nii','.nrrd','.nrrd.gz','.gipl','.gipl.gz']
    file_list = GetListFiles(folder_path, file_extension)

    patients = {}
    
    for file in file_list:
        basename = os.path.basename(file)
        patient = basename.split('_Or')[0].split('_MAND')[0].split('_MAX')[0].split('_CB')[0].split('_T1')[0].split('_T2')[0].split('.')[0]
        
        if patient not in patients:
            patients[patient] = {}
        
        if 'Seg' in basename or 'seg' in basename:
            patients[patient]['seg'+time_point] = file
        
        else:
            patients[patient]['scan'+time_point] = file
    
    return patients

def MergeDicts(dict1,dict2):
    """Merge t1 and t2 dictionaries for each patient"""
    patients = {}
    for patient in dict1:
        patients[patient] = dict1[patient]
        try:
            patients[patient].update(dict2[patient])
        except KeyError:
            continue
    return patients

def search(path,*args):
    """
    Return a dictionary with args element as key and a list of file in path directory finishing by args extension for each key

    Example:
    args = ('json',['.nii.gz','.nrrd'])
    return:
        {
            'json' : ['path/a.json', 'path/b.json','path/c.json'],
            '.nii.gz' : ['path/a.nii.gz', 'path/b.nii.gz']
            '.nrrd.gz' : ['path/c.nrrd']
        }
    """
    arguments=[]
    for arg in args:
        if type(arg) == list:
            arguments.extend(arg)
        else:
            arguments.append(arg)
    return {key: sorted([i for i in iglob(os.path.normpath("/".join([path,'**','*'])),recursive=True) if i.endswith(key)]) for key in arguments}

def ResampleImage(image, transform):
    '''
    Resample image using SimpleITK
    
    Parameters
    ----------
    image : SimpleITK.Image
        Image to be resampled
    target : SimpleITK.Image
        Target image
    transform : SimpleITK transform
        Transform to be applied to the image.
        
    Returns
    -------
    SimpleITK image
        Resampled image.
    '''
    # Create resampler
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(transform)

    # Resample image
    resampled_image = resampler.Execute(image)

    return resampled_image

def SitkRegTransformNonGrowing(fixed_image, moving_image,outpath):
    """Get the transform to register two images for a non growing patient using SimpleITK registration method"""

    # Register images
    initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    registration_method = sitk.ImageRegistrationMethod()
    sitk.WriteTransform(initial_transform, outpath+'/initial_transform.tfm')
    registration_method.SetInitialTransform(initial_transform)
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    # registration_method.SetMetricAsMeanSquares()
    # registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    # registration_method.SetMetricSamplingPercentage(0.01)
    registration_method.SetInterpolator(sitk.sitkLinear)

    registration_method.SetOptimizerAsGradientDescent(learningRate=1, numberOfIterations=1000, estimateLearningRate=registration_method.EachIteration)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])

    final_transform_1 = registration_method.Execute(fixed_image, moving_image)
    sitk.WriteTransform(final_transform_1, outpath+'/final_transform_1.tfm')

    # Fine Tuning
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetInitialTransform(final_transform_1)
    registration_method.SetOptimizerAsGradientDescent(learningRate=1, numberOfIterations=100)
    registration_method.SetOptimizerScalesFromPhysicalShift()
    registration_method.SetMetricAsMeanSquares()
    registration_method.SetInterpolator(sitk.sitkLinear)

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[0])

    final_transform_2 = registration_method.Execute(fixed_image, moving_image)
    sitk.WriteTransform(final_transform_2, outpath+'/final_transform_2.tfm')
    return final_transform_1

def SitkRegTransformGrowing(fixed_image, moving_image):
    """Get the transform to register two images for a growing patient using SimpleITK registration method"""

    # Initial Alignment
    initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)

    # Non Linear Registration
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMeanSquares()
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=1e-4, numberOfIterations=100000)#, maximumNumberOfCorrections=5, maximumNumberOfFunctionEvaluations=1000, costFunctionConvergenceFactor=1e+7)
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])

    final_transform = registration_method.Execute(sitk.Cast(fixed_image,sitk.sitkFloat32), sitk.Cast(moving_image,sitk.sitkFloat32))

    # Fine Tuning
    registration_method.SetOptimizerAsGradientDescent(learningRate=1e-4, numberOfIterations=100000)
    registration_method.SetInitialTransform(final_transform)

    final_transform = registration_method.Execute(sitk.Cast(fixed_image,sitk.sitkFloat32), sitk.Cast(moving_image,sitk.sitkFloat32))

    return final_transform

def applyMask(image, mask):
    """Apply a mask to an image."""
    # Cast the image to float32
    image = sitk.Cast(image, sitk.sitkFloat32)

    return sitk.Mask(image, mask)

def VoxelBasedRegistration(fixed_image_path,moving_image_path,fixed_seg_path,moving_seg_path,output_dir,patient):
    outpath = os.path.join(output_dir,patient+'_OutReg')

    if os.path.exists(outpath):
        shutil.rmtree(outpath)

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Copy T1 and T2 images to output directory
    shutil.copyfile(fixed_image_path, os.path.join(outpath,patient+'_T1.nii.gz'))
    shutil.copyfile(moving_image_path, os.path.join(outpath,patient+'_T2.nii.gz'))
    
    # Read images and segmentations
    fixed_image = sitk.ReadImage(fixed_image_path)
    fixed_seg = sitk.ReadImage(fixed_seg_path)
    moving_image = sitk.ReadImage(moving_image_path)
    moving_seg = sitk.ReadImage(moving_seg_path)

    # Print sizes
    print('Fixed image size: ', fixed_image.GetSize())
    print('Moving image size: ', moving_image.GetSize())

    print('Fixed image origin: ', fixed_image.GetOrigin())

    # Apply mask to images
    fixed_image_masked = applyMask(fixed_image, fixed_seg)
    moving_image_masked = applyMask(moving_image, moving_seg)
    
    # Register images
    # tic = time.time()
    transform = SitkRegTransformNonGrowing(fixed_image_masked, moving_image_masked,outpath)
    # print('Registration time: ', round(time.time() - tic,2),'s')
    # sitk.WriteTransform(transform, os.path.join(outpath,patient+'_matrix.tfm'))

    # Resample images and segmentations using the final transform
    # tic = time.time()
    resample_t2 = ResampleImage(moving_image, transform)
    print('Resampled image size: ', resample_t2.GetSize())
    print('Resampled image origin: ', resample_t2.GetOrigin())

    resample_t2_seg = ResampleImage(moving_seg, transform)
    # sitk.WriteImage(resample_t2, os.path.join(outpath,patient+'_ScanReg.nii.gz'))
    # sitk.WriteImage(resample_t2_seg, os.path.join(outpath,patient+'_SegReg.nii.gz'))
    # print('Resampling time: ', round(time.time() - tic,2),'s')

def main(args):
    t1_folder,t2_folder, output_dir = args.t1_folder, args.t2_folder, args.output_dir

    t1_patients = GetPatients(t1_folder, time_point='T1')
    t2_patients = GetPatients(t2_folder, time_point='T2')

    # merge dictionaries for each patient
    patients = MergeDicts(t1_patients,t2_patients)

    for patient,data in patients.items():
        print("="*70)
        print("Working on patient {}".format(patient))
        try:
            VoxelBasedRegistration(data['scanT1'],data['scanT2'],data['segT1'],data['segT2'],output_dir,patient)
        except KeyError:
            print("Patient {} does not have both T1 and T2 scans".format(patient))
            continue

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Voxel Based Registration')
    parser.add_argument('--t1_folder', type=str, help='Path to folder containing input T1 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/Selene/T1')
    parser.add_argument('--t2_folder', type=str, help='Path to folder containing input T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/Selene/T2_not_registeredOr')
    parser.add_argument('--output_dir', type=str, help='Path to folder containing output register T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/Selene/T2_auto_reg')

    args = parser.parse_args()

    main(args)