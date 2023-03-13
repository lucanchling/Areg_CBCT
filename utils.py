"""
8888888 888b     d888 8888888b.   .d88888b.  8888888b.  88888888888  .d8888b.  
  888   8888b   d8888 888   Y88b d88P" "Y88b 888   Y88b     888     d88P  Y88b 
  888   88888b.d88888 888    888 888     888 888    888     888     Y88b.      
  888   888Y88888P888 888   d88P 888     888 888   d88P     888      "Y888b.   
  888   888 Y888P 888 8888888P"  888     888 8888888P"      888         "Y88b. 
  888   888  Y8P  888 888        888     888 888 T88b       888           "888 
  888   888   "   888 888        Y88b. .d88P 888  T88b      888     Y88b  d88P 
8888888 888       888 888         "Y88888P"  888   T88b     888      "Y8888P" 
"""
import SimpleITK as sitk
import numpy as np
import time,shutil
from glob import iglob
import os
import argparse
from matplotlib import pyplot as plt
import json
from tqdm import tqdm

"""
8888888888 8888888 888      8888888888  .d8888b.  
888          888   888      888        d88P  Y88b 
888          888   888      888        Y88b.      
8888888      888   888      8888888     "Y888b.   
888          888   888      888            "Y88b. 
888          888   888      888              "888 
888          888   888      888        Y88b  d88P 
888        8888888 88888888 8888888888  "Y8888P"
"""

def GetListFiles(folder_path, file_extension):
    """Return a list of files in folder_path finishing by file_extension"""
    file_list = []
    for extension_type in file_extension:
        file_list += search(folder_path,file_extension)[extension_type]
    return file_list

def GetPatients(folder_path, time_point='T1', segmentationType='CB'):
    """Return a dictionary with patient id as key"""

    file_extension = ['.nii.gz','.nii','.nrrd','.nrrd.gz','.gipl','.gipl.gz']
    json_extension = ['.json']
    file_list = GetListFiles(folder_path, file_extension+json_extension)

    patients = {}
    
    for file in file_list:
        basename = os.path.basename(file)
        patient = basename.split('_Or')[0].split('_ScanReg')[0].split('_OR')[0].split('_MAND')[0].split('_MAX')[0].split('_CB')[0].split('_T1')[0].split('_T2')[0].split('.')[0].split('_lm')[0]
        
        if patient not in patients:
            patients[patient] = {}
        
        if True in [i in basename for i in file_extension]:
            if segmentationType+'MASK' in basename:
                patients[patient]['seg'+time_point] = file
            
            if 'MASK' not in basename:
                patients[patient]['scan'+time_point] = file

        else:
            if time_point == 'T2':
                patients[patient]['lm'+time_point] = file

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

"""
888             d8888 888b    888 8888888b.  888b     d888        d8888 8888888b.  888    d8P  
888            d88888 8888b   888 888  "Y88b 8888b   d8888       d88888 888   Y88b 888   d8P   
888           d88P888 88888b  888 888    888 88888b.d88888      d88P888 888    888 888  d8P    
888          d88P 888 888Y88b 888 888    888 888Y88888P888     d88P 888 888   d88P 888d88K     
888         d88P  888 888 Y88b888 888    888 888 Y888P 888    d88P  888 8888888P"  8888888b    
888        d88P   888 888  Y88888 888    888 888  Y8P  888   d88P   888 888 T88b   888  Y88b   
888       d8888888888 888   Y8888 888  .d88P 888   "   888  d8888888888 888  T88b  888   Y88b  
88888888 d88P     888 888    Y888 8888888P"  888       888 d88P     888 888   T88b 888    Y88b 
"""

def applyTransformLandmarks(landmarks, transform):
    """Apply a transform to a set of landmarks."""
    copy = landmarks.copy()
    for lm, pt in landmarks.items():
        copy[lm] = transform.TransformPoint(pt)
    return copy

def GenControlePoint(landmarks):
    lm_lst = []
    false = False
    true = True
    id = 0
    for landmark,data in landmarks.items():
        id+=1
        controle_point = {
            "id": str(id),
            "label": landmark,
            "description": "",
            "associatedNodeID": "",
            "position": [data[0], data[1], data[2]],
            "orientation": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            "selected": true,
            "locked": true,
            "visibility": true,
            "positionStatus": "defined"
        }
        lm_lst.append(controle_point)

    return lm_lst

def WriteJson(landmarks,out_path):
    false = False
    true = True
    file = {
    "@schema": "https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/Markups/Resources/Schema/markups-schema-v1.0.0.json#",
    "markups": [
        {
            "type": "Fiducial",
            "coordinateSystem": "LPS",
            "locked": false,
            "labelFormat": "%N-%d",
            "controlPoints": GenControlePoint(landmarks),
            "measurements": [],
            "display": {
                "visibility": false,
                "opacity": 1.0,
                "color": [0.4, 1.0, 0.0],
                "color": [0.5, 0.5, 0.5],
                "selectedColor": [0.26666666666666669, 0.6745098039215687, 0.39215686274509806],
                "propertiesLabelVisibility": false,
                "pointLabelsVisibility": true,
                "textScale": 2.0,
                "glyphType": "Sphere3D",
                "glyphScale": 2.0,
                "glyphSize": 5.0,
                "useGlyphScale": true,
                "sliceProjection": false,
                "sliceProjectionUseFiducialColor": true,
                "sliceProjectionOutlinedBehindSlicePlane": false,
                "sliceProjectionColor": [1.0, 1.0, 1.0],
                "sliceProjectionOpacity": 0.6,
                "lineThickness": 0.2,
                "lineColorFadingStart": 1.0,
                "lineColorFadingEnd": 10.0,
                "lineColorFadingSaturation": 1.0,
                "lineColorFadingHueOffset": 0.0,
                "handlesInteractive": false,
                "snapMode": "toVisibleSurface"
            }
        }
    ]
    }
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(file, f, ensure_ascii=False, indent=4)

    f.close

def LoadOnlyLandmarks(ldmk_path, ldmk_list=None):
    """
    Load landmarks from json file without using the img as input
    
    Parameters
    ----------
    ldmk_path : str
        Path to the json file
    gold : bool, optional
        If True, load gold standard landmarks, by default False
    
    Returns
    -------
    dict
        Dictionary of landmarks
    
    Raises
    ------
    ValueError
        If the json file is not valid
    """
    with open(ldmk_path) as f:
        data = json.load(f)
    
    markups = data["markups"][0]["controlPoints"]
    
    landmarks = {}
    for markup in markups:
        try:
            lm_ph_coord = np.array([markup["position"][0],markup["position"][1],markup["position"][2]])
            #lm_coord = ((lm_ph_coord - origin) / spacing).astype(np.float16)
            lm_coord = lm_ph_coord.astype(np.float64)
            landmarks[markup["label"]] = lm_coord
        except:
            continue
    if ldmk_list is not None:
        return {key:landmarks[key] for key in ldmk_list if key in landmarks.keys()}
    
    return landmarks

"""
8888888 888b     d888        d8888  .d8888b.  8888888888  .d8888b.  
  888   8888b   d8888       d88888 d88P  Y88b 888        d88P  Y88b 
  888   88888b.d88888      d88P888 888    888 888        Y88b.      
  888   888Y88888P888     d88P 888 888        8888888     "Y888b.   
  888   888 Y888P 888    d88P  888 888  88888 888            "Y88b. 
  888   888  Y8P  888   d88P   888 888    888 888              "888 
  888   888   "   888  d8888888888 Y88b  d88P 888        Y88b  d88P 
8888888 888       888 d88P     888  "Y8888P88 8888888888  "Y8888P" 
"""

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

def applyMask(image, mask):
    """Apply a mask to an image."""
    # Cast the image to float32
    # image = sitk.Cast(image, sitk.sitkFloat32)

    return sitk.Mask(image, mask)


"""
8888888b.  8888888888  .d8888b.  8888888  .d8888b.  88888888888 
888   Y88b 888        d88P  Y88b   888   d88P  Y88b     888     
888    888 888        888    888   888   Y88b.          888     
888   d88P 8888888    888          888    "Y888b.       888     
8888888P"  888        888  88888   888       "Y88b.     888     
888 T88b   888        888    888   888         "888     888     
888  T88b  888        Y88b  d88P   888   Y88b  d88P     888     
888   T88b 8888888888  "Y8888P88 8888888  "Y8888P"      888
"""

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

def MatrixRetrieval(TransformParameterMap):
    Transforms = []
    Matrix = np.zeros((4,4))
    for ParameterMap in TransformParameterMap:
        if ParameterMap['Transform'][0] == 'AffineTransform':
            matrix = [float(i) for i in ParameterMap['TransformParameters']]
            # Convert to a sitk transform
            transform = sitk.AffineTransform(3)
            transform.SetParameters(matrix)
            Transforms.append(transform)

        elif ParameterMap['Transform'][0] == 'EulerTransform':
            A = [float(i) for i in ParameterMap['TransformParameters'][0:3]]
            B = [float(i) for i in ParameterMap['TransformParameters'][3:6]]
            # Convert to a sitk transform
            transform = sitk.Euler3DTransform()
            transform.SetRotation(angleX=A[0], angleY=A[1], angleZ=A[2])
            transform.SetTranslation(B)
            Transforms.append(transform)
    
    # Create a composite transform
    final_transform = sitk.Transform()
    for transform in Transforms:
        final_transform.AddTransform(transform)
    
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
    elastixImageFilter.SetParameter("MaximumNumberOfIterations", "100")
    # elastixImageFilter.SetParameter("NumberOfResolutions", "4")
    
    tic = time.time()
    elastixImageFilter.Execute()
    # print('Process Time: {} sec'.format(round(time.time()-tic,2)))

    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    return resultImage, transformParameterMap

def VoxelBasedRegistration(fixed_image_path,moving_image_path,fixed_seg_path,moving_seg_path,outpath,patient):
    if os.path.exists(outpath):
        shutil.rmtree(outpath)

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Copy T1 and T2 images to output directory
    # shutil.copyfile(fixed_image_path, os.path.join(outpath,patient+'_T1.nii.gz'))
    # shutil.copyfile(moving_image_path, os.path.join(outpath,patient+'_T2.nii.gz'))
    
    # Read images and segmentations
    fixed_image = sitk.ReadImage(fixed_image_path)
    fixed_seg = sitk.ReadImage(fixed_seg_path)
    moving_image = sitk.ReadImage(moving_image_path)
    moving_seg = sitk.ReadImage(moving_seg_path)

    # Apply mask to images
    fixed_image_masked = applyMask(fixed_image, fixed_seg)
    # moving_image_masked = applyMask(moving_image, moving_seg)

    # Register images
    # tic = time.time()
    resample_t2, TransformParamMap = SimpleElastixReg(fixed_image_masked, moving_image,outpath)
    transform = MatrixRetrieval(TransformParamMap)
    # print('Registration time: ', round(time.time() - tic,2),'s')
    sitk.WriteTransform(transform, os.path.join(outpath,patient+'_matrix.tfm'))

    # Resample images and segmentations using the final transform
    # tic = time.time()
    resample_t2 = sitk.Cast(ResampleImage(moving_image, transform),sitk.sitkInt16)
    resample_t2_seg = ResampleImage(moving_seg, transform)

    
    # Compare segmentations
    # DiceStart,JaccardStart,HausdorfStart = CompareScans(fixed_seg, moving_seg)
    # DiceEnd,JaccardEnd,HausdorfEnd = CompareScans(fixed_seg, resample_t2_seg)
    # print("Delta Dice {} | Delta Jaccard {} | Delta Hausdorf {}".format(round(DiceEnd-DiceStart,5),round(JaccardEnd-JaccardStart,5),round(HausdorfEnd-HausdorfStart,5)))
    # sitk.WriteImage(sitk.Cast(resample_t2,sitk.sitkInt16), os.path.join(outpath,patient+'_ScanReg.nii.gz'))
    # print('Resampling time: ', round(time.time() - tic,2),'s')

    return transform, resample_t2, resample_t2_seg


"""
888     888 88888888888 8888888 888       .d8888b.  
888     888     888       888   888      d88P  Y88b 
888     888     888       888   888      Y88b.      
888     888     888       888   888       "Y888b.   
888     888     888       888   888          "Y88b. 
888     888     888       888   888            "888 
Y88b. .d88P     888       888   888      Y88b  d88P 
 "Y88888P"      888     8888888 88888888  "Y8888P"
"""

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

def translate(shortname):
    dic = {'CB': 'Cranial Base', 'MAND': 'Mandible', 'MAX': 'Maxilla'}
    return dic[shortname]
