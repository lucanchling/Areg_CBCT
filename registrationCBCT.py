import SimpleITK as sitk
import numpy as np
import time,shutil
from glob import iglob
import os
import argparse
from matplotlib import pyplot as plt


def GetListFiles(folder_path, file_extension):
    """Return a list of files in folder_path finishing by file_extension"""
    file_list = []
    for extension_type in file_extension:
        file_list += search(folder_path,file_extension)[extension_type]
    return file_list

def GetPatients(folder_path, time_point='T1', segmentationType='CB'):
    """Return a dictionary with patient id as key"""

    file_extension = ['.nii.gz','.nii','.nrrd','.nrrd.gz','.gipl','.gipl.gz']
    file_list = GetListFiles(folder_path, file_extension)

    patients = {}
    
    for file in file_list:
        basename = os.path.basename(file)
        patient = basename.split('_Or')[0].split('_MAND')[0].split('_MAX')[0].split('_CB')[0].split('_T1')[0].split('_T2')[0].split('.')[0]
        
        if patient not in patients:
            patients[patient] = {}
        
        if segmentationType+'MASK' in basename:
            patients[patient]['seg'+time_point] = file
        
        if 'MASK' not in basename:
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

def CorrectHisto(input_img,min_porcent=0.01,max_porcent = 0.99, i_min=-1500, i_max=4000):
    input_img = sitk.Cast(input_img, sitk.sitkFloat32)
    img = sitk.GetArrayFromImage(input_img)


    img_min = np.min(img)
    img_max = np.max(img)
    img_range = img_max - img_min
    # print(img_min,img_max,img_range)

    definition = 1000
    histo = np.histogram(img,definition)
    cum = np.cumsum(histo[0])
    cum = cum - np.min(cum)
    cum = cum / np.max(cum)

    res_high = list(map(lambda i: i> max_porcent, cum)).index(True)
    res_max = (res_high * img_range)/definition + img_min

    res_low = list(map(lambda i: i> min_porcent, cum)).index(True)
    res_min = (res_low * img_range)/definition + img_min

    res_min = max(res_min,i_min)
    res_max = min(res_max,i_max)


    # print(res_min,res_min)

    img = np.where(img > res_max, res_max,img)
    img = np.where(img < res_min, res_min,img)

    image = sitk.GetImageFromArray(img)
    image.CopyInformation(input_img)

    return image
    
def demons_registration(fixed_image, moving_image):
    registration_method = sitk.ImageRegistrationMethod()

    # Create initial identity transformation.
    transform_to_displacement_field_filter = sitk.TransformToDisplacementFieldFilter()
    transform_to_displacement_field_filter.SetReferenceImage(fixed_image)
    # The image returned from the initial_transform_filter is transferred to the transform and cleared out.
    initial_transform = sitk.DisplacementFieldTransform(transform_to_displacement_field_filter.Execute(sitk.Transform()))
    
    # Regularization (update field - viscous, total field - elastic).
    initial_transform.SetSmoothingGaussianOnUpdate(varianceForUpdateField=0.0, varianceForTotalField=2.0) 
    
    registration_method.SetInitialTransform(initial_transform)

    registration_method.SetMetricAsDemons(10) #intensities are equal if the difference is less than 10HU
        
    # Multi-resolution framework.            
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = [8,4,0])    

    registration_method.SetInterpolator(sitk.sitkLinear)
    # If you have time, run this code using the ConjugateGradientLineSearch, otherwise run as is.   
    # registration_method.SetOptimizerAsConjugateGradientLineSearch(learningRate=1.0, numberOfIterations=20, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=1)#, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()
        
    return registration_method.Execute(fixed_image, moving_image)

def MeanSquaresReg(fixed_image, moving_image, lr=1.0, nbIterations=500, initial_transform=None, sampling=None,shrinkFactors = [4,2,1],smoothingSigmas = [8,4,0]):
    
    registration_method = sitk.ImageRegistrationMethod()
    
    if initial_transform is None:
        initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    registration_method.SetInitialTransform(initial_transform)

    registration_method.SetMetricAsMeanSquares()
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=lr, numberOfIterations=nbIterations)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    if sampling is not None:
        registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
        registration_method.SetMetricSamplingPercentage(sampling)
    else:
        registration_method.SetMetricSamplingStrategy(registration_method.NONE)

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = shrinkFactors)
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = smoothingSigmas)

    return registration_method.Execute(fixed_image, moving_image)

def CorrelationReg(fixed_image, moving_image, lr=1.0, nbIterations=500, initial_transform=None, sampling=None,shrinkFactors = [4,2,1],smoothingSigmas = [8,4,0]):
    
    registration_method = sitk.ImageRegistrationMethod()
    
    if initial_transform is None:
        initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    registration_method.SetInitialTransform(initial_transform)

    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=lr, numberOfIterations=nbIterations)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    if sampling is not None:
        registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
        registration_method.SetMetricSamplingPercentage(sampling)
    else:
        registration_method.SetMetricSamplingStrategy(registration_method.NONE)

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = shrinkFactors)
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = smoothingSigmas)

    return registration_method.Execute(fixed_image, moving_image)

def MattesMutualInformationReg(fixed_image, moving_image,numberOfHistogramBins=50, lr=1.0, nbIterations=500, initial_transform=None, sampling=None,shrinkFactors = [4,2,1],smoothingSigmas = [8,4,0]):
    
    registration_method = sitk.ImageRegistrationMethod()
    
    if initial_transform is None:
        initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    registration_method.SetInitialTransform(initial_transform)

    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=numberOfHistogramBins)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=lr, numberOfIterations=nbIterations)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    if sampling is not None:
        registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
        registration_method.SetMetricSamplingPercentage(sampling)
    else:
        registration_method.SetMetricSamplingStrategy(registration_method.NONE)

    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = shrinkFactors)
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = smoothingSigmas)

    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    return registration_method.Execute(fixed_image, moving_image)

def SitkRegTransformNonGrowing(fixed_image, moving_image,outpath):
    """Get the transform to register two images for a non growing patient using SimpleITK registration method"""

    # First Registration Step
    tic_global = time.time()
    tic = time.time()
    final_transform = MattesMutualInformationReg(fixed_image, moving_image, sampling=0.01, nbIterations=10000, lr=1, shrinkFactors=[4,2,1],smoothingSigmas=[2,1,0])
    # sitk.WriteTransform(final_transform, outpath+'/trans1.tfm')
    # sitk.WriteTransform(final_transform_1, outpath+'/final_transform_1.tfm')
    # print('First Registration: {} sec'.format(round(time.time()-tic,2)))

    # Second Registration Step
    tic = time.time()
    final_transform = MattesMutualInformationReg(fixed_image, moving_image, sampling=0.01, nbIterations=20000, lr=1e-4, shrinkFactors=[1], smoothingSigmas=[0], initial_transform=final_transform)
    # sitk.WriteTransform(final_transform, outpath+'/trans2.tfm')
    # print('Second Registration: {} sec'.format(round(time.time()-tic,2)))
    
    # Third Registration Step
    tic = time.time()
    # final_transform = CorrelationReg(fixed_image, moving_image, sampling=0.01, nbIterations=100, lr=1e-4, shrinkFactors=[1], smoothingSigmas=[0], initial_transform=final_transform)
    # print('Third Registration: {} sec'.format(round(time.time()-tic,2)))
    
    # # Fourth Registration Step
    # tic = time.time()
    # final_transform = MattesMutualInformationReg(fixed_image, moving_image, sampling=0.01, nbIterations=100, lr=1e-2, shrinkFactors=[1], smoothingSigmas=[0], initial_transform=final_transform)
    # print('Fourth Registration: {} sec'.format(round(time.time()-tic,2)))

    # # Fifth Registration Step
    # tic = time.time()
    # final_transform = CorrelationReg(fixed_image, moving_image, sampling=0.01, nbIterations=25, lr=1e-3, shrinkFactors=[1], smoothingSigmas=[0], initial_transform=final_transform)
    # print('Fifth Registration: {} sec'.format(round(time.time()-tic,2)))

    tic = time.time()
    # final_transform = CorrelationReg(fixed_image, moving_image, sampling=0.1, nbIterations=100, lr=1e-2, shrinkFactors=[1], smoothingSigmas=[0], initial_transform=final_transform)
    
    sitk.WriteTransform(final_transform, outpath+'/'+os.path.basename(outpath).split('_OutReg')[0]+'_trans.tfm')
    print('Total Process Time: {} sec'.format(round(time.time()-tic_global,2)))
    
    return final_transform

def SimpleElastixReg(fixed_image, moving_image, outpath):
    """Get the transform to register two images using SimpleElastix registration method"""
    
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.LogToConsoleOff()
    elastixImageFilter.SetFixedImage(fixed_image)
    elastixImageFilter.SetMovingImage(moving_image)

    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("affine"))
    parameterMapVector.append(sitk.GetDefaultParameterMap("rigid"))
    elastixImageFilter.SetParameterMap(parameterMapVector)
    
    elastixImageFilter.SetParameter("ErodeMask", "true")
    # elastixImageFilter.SetParameter("MaximumNumberOfIterations", "5000")
    # elastixImageFilter.SetParameter("NumberOfResolutions", "4")
    
    tic = time.time()
    elastixImageFilter.Execute()
    print('Process Time: {} sec'.format(round(time.time()-tic,2)))

    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    return resultImage

def printHist(fixed_image, moving_image, fixed_masked_image, moving_masked_image):
    """Print the histograms of the fixed and moving images"""
    
    fixed_image_array = sitk.GetArrayFromImage(fixed_image)
    moving_image_array = sitk.GetArrayFromImage(moving_image)
    fixed_masked_image_array = sitk.GetArrayFromImage(fixed_masked_image)
    moving_masked_image_array = sitk.GetArrayFromImage(moving_masked_image)
    
    
    plt.figure()
    plt.hist(fixed_image_array.flatten(), bins=50, alpha=0.5, label='fixed')
    plt.hist(moving_image_array.flatten(), bins=50, alpha=0.5, label='moving')
    plt.hist(fixed_masked_image_array.flatten(), bins=50, alpha=0.5, label='fixed_masked')
    plt.hist(moving_masked_image_array.flatten(), bins=50, alpha=0.5, label='moving_masked')
    plt.legend(loc='upper right')
    plt.show()

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
    # image = sitk.Cast(image, sitk.sitkFloat32)

    return sitk.Mask(image, mask)

def VoxelBasedRegistration(fixed_image_path,moving_image_path,fixed_seg_path,moving_seg_path,output_dir,patient):
    outpath = os.path.join(output_dir,patient+'_OutReg')

    if os.path.exists(outpath):
        shutil.rmtree(outpath)

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Copy T1 and T2 images to output directory
    shutil.copyfile(fixed_image_path, os.path.join(outpath,patient+'_T1.nii.gz'))
    # shutil.copyfile(moving_image_path, os.path.join(outpath,patient+'_T2.nii.gz'))
    
    # Read images and segmentations
    fixed_image = sitk.ReadImage(fixed_image_path)
    fixed_seg = sitk.ReadImage(fixed_seg_path)
    moving_image = sitk.ReadImage(moving_image_path)
    moving_seg = sitk.ReadImage(moving_seg_path)

    # Print Intensity Range
    # print('Fixed image intensity range: ', np.min(array_fixed_image), np.max(array_fixed_image))
    # print('Moving image intensity range: ', np.min(array_mov_image), np.max(array_mov_image))
    # Print sizes
    # print('Fixed image size: ', fixed_image.GetSize())
    # print('Moving image size: ', moving_image.GetSize())
    # print('Fixed image spacing: ', fixed_image.GetSpacing())
    # print('Moving image spacing: ', moving_image.GetSpacing())

    # print('Fixed image origin: ', fixed_image.GetOrigin())

    # Apply mask to images
    fixed_image_masked = applyMask(fixed_image, fixed_seg)
    moving_image_masked = applyMask(moving_image, moving_seg)

    # print('Fixed image masked size: ', fixed_image_masked.GetSize())
    # print('Moving image masked size: ', moving_image_masked.GetSize())
    # print('Fixed image masked spacing: ', fixed_image_masked.GetSpacing())
    # print('Moving image masked spacing: ', moving_image_masked.GetSpacing())

    # printHist(fixed_image, moving_image, fixed_image_masked, moving_image_masked)

    # Print Intensity Range
    # print('Fixed image masked intensity range: ', np.min(array_fixed_image_masked), np.max(array_fixed_image_masked))
    # print('Moving image masked intensity range: ', np.min(array_mov_image_masked), np.max(array_mov_image_masked))

    # printHist(fixed_image_masked, moving_image_masked)
    
    # Register images
    # tic = time.time()
    resample_t2 = SimpleElastixReg(fixed_image_masked, moving_image,outpath)
    # print('Registration time: ', round(time.time() - tic,2),'s')
    # sitk.WriteTransform(transform, os.path.join(outpath,patient+'_matrix.tfm'))

    # Resample images and segmentations using the final transform
    # tic = time.time()
    # resample_t2 = ResampleImage(moving_image_masked, transform)
    # print('Resampled image size: ', resample_t2.GetSize())
    # print('Resampled image origin: ', resample_t2.GetOrigin())

    # resample_t2_seg = ResampleImage(moving_seg, transform)
    
    # Compare segmentations
    # DiceStart,JaccardStart,HausdorfStart = CompareScans(fixed_seg, moving_seg)
    # DiceEnd,JaccardEnd,HausdorfEnd = CompareScans(fixed_seg, resample_t2_seg)
    # print("Delta Dice {} | Delta Jaccard {} | Delta Hausdorf {}".format(round(DiceEnd-DiceStart,5),round(JaccardEnd-JaccardStart,5),round(HausdorfEnd-HausdorfStart,5)))
    sitk.WriteImage(sitk.Cast(resample_t2,sitk.sitkInt16), os.path.join(outpath,patient+'_ScanReg.nii.gz'))
    # sitk.WriteImage(resample_t2_seg, os.path.join(outpath,patient+'_SegReg.nii.gz'))
    # print('Resampling time: ', round(time.time() - tic,2),'s')

def CompareScans(im1, im2):
    """Compare two segmentations with Dice similarity coefficient, Jaccard similarity coefficient, Intersection over Union, and Hausdorff distance."""
    # Cast the images to float32
    im1 = sitk.Cast(im1, sitk.sitkInt16)
    im2 = sitk.Cast(im2, sitk.sitkInt16)

    # Compute Dice similarity coefficient
    OverlaFilter = sitk.LabelOverlapMeasuresImageFilter()
    OverlaFilter.Execute(im1, im2)
    dice = OverlaFilter.GetDiceCoefficient()        # Goal: 0.7
    jaccard = OverlaFilter.GetJaccardCoefficient()  # Goal: 0.6

    # Compute Hausdorff distance
    hausdorff = sitk.HausdorffDistanceImageFilter()
    hausdorff.Execute(im1, im2)
    hausdorff = hausdorff.GetHausdorffDistance()
    return dice,jaccard,hausdorff

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
            # print("T1 scan path: {} \nT2 scan path: {}".format(os.path.basename(data['scanT1']),os.path.basename(data['scanT2'])))
            # print("T1 seg path: {} \nT2 seg path: {}".format(os.path.basename(data['segT1']),os.path.basename(data['segT2'])))
            VoxelBasedRegistration(data['scanT1'],data['scanT2'],data['segT1'],data['segT2'],output_dir,patient)
        except KeyError:
            print("Patient {} does not have both T1 and T2 scans".format(patient))
            continue

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Voxel Based Registration')
    parser.add_argument('--t1_folder', type=str, help='Path to folder containing input T1 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/T1/')
    parser.add_argument('--t2_folder', type=str, help='Path to folder containing input T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/T2/')
    parser.add_argument('--output_dir', type=str, help='Path to folder containing output register T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/OUT/')

    args = parser.parse_args()

    main(args)