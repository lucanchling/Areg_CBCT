import SimpleITK as sitk
import itk
from MergeTransform import ComputeFinalMatrix
from utils import VoxelBasedRegistration

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


def MatrixRetrieval(TransformParameterMapObject):
    """Retrieve the matrix from the transform parameter map"""
    Transforms = []
    ParameterMap = TransformParameterMapObject.GetParameterMap(0)

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

    
    return Transforms

def ElastixApprox(fixed_image, moving_image):
    elastix_object = itk.ElastixRegistrationMethod.New(fixed_image, moving_image)
    
    # ParameterMap 
    parameter_object = itk.ParameterObject.New()
    default_rigid_parameter_map = parameter_object.GetDefaultParameterMap('rigid')
    parameter_object.AddParameterMap(default_rigid_parameter_map)
    parameter_object.SetParameter("WriteResultImage", "false")
    
    elastix_object.SetParameterObject(parameter_object)

    # Additional parameters
    elastix_object.SetLogToConsole(False)

    # Execute registration
    elastix_object.UpdateLargestPossibleRegion()
    
    TransParamObj = elastix_object.GetTransformParameterObject()
    
    return TransParamObj

def ElastixReg(fixed_image, moving_image, fixed_mask, initial_transform=None):

    elastix_object = itk.ElastixRegistrationMethod.New(fixed_image, moving_image)
    # elastix_object.SetFixedMask(fixed_mask)

    # ParameterMap 
    parameter_object = itk.ParameterObject.New()
    default_rigid_parameter_map = parameter_object.GetDefaultParameterMap('rigid')
    parameter_object.AddParameterMap(default_rigid_parameter_map)
    parameter_object.SetParameter("ErodeMask", "true")
    parameter_object.SetParameter("WriteResultImage", "false")
    parameter_object.SetParameter("MaximumNumberOfIterations", "10000")
    parameter_object.SetParameter("NumberOfResolutions", "1")
    parameter_object.SetParameter("NumberOfSpatialSamples", "10000")
    # parameter_object.SetParameter("MaximumNumberOfSamplingAttempts", "25")

    elastix_object.SetParameterObject(parameter_object)
    if initial_transform is not None:
        elastix_object.SetInitialTransformParameterObject(initial_transform)

    # Additional parameters
    elastix_object.SetLogToConsole(False)

    # Execute registration
    elastix_object.UpdateLargestPossibleRegion()
    
    TransParamObj = elastix_object.GetTransformParameterObject()
    
    return TransParamObj



if __name__ == "__main__":
    fixed_path = '/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/JJ/Approach2/CB/T1/Pat1/Pat1_Cl_II_ScanT1_Or.gipl.gz'
    fixed_mask_path = '/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/JJ/Approach2/CB/T1/Pat1/Pat1_Cl_II_SegmT1_Or.gipl.gz'
    moving_path = '/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/JJ/Approach3/T2_Center/Pat1/Pat1_Cl_II_ScanT2.gipl.gz'
    # transformSITK = VoxelBasedRegistration(fixed_path,moving_path,fixed_mask_path,approx=True)
    # sitk.WriteTransform(transformSITK,'/home/luciacev/Downloads/TransformSITK.tfm')
    fixed_image = itk.imread(fixed_path,itk.F)
    moving_image = itk.imread(moving_path,itk.F)

    fixed_mask = itk.imread(fixed_mask_path,itk.UC)

    tmp_image = sitk.ReadImage(fixed_path)
    sitk_masked_image = sitk.Mask(tmp_image,sitk.ReadImage(fixed_mask_path))
    # Copy the information from the original image
    sitk_masked_image.CopyInformation(tmp_image)

    sitk.WriteImage(sitk_masked_image,'/home/luciacev/Downloads/masked_image.nii.gz')
    
    itk_masked_image = itk.imread('/home/luciacev/Downloads/masked_image.nii.gz',itk.F)
    # Mask the image
    # mask_filter = itk.MaskImageFilter.New(fixed_image)
    # mask_filter.SetMaskImage(fixed_mask)
    # mask_filter.SetOutsideValue(0)
    # fixed_image_masked = mask_filter.GetOutput()


    inittransform = ElastixApprox(fixed_image,moving_image)
    init_transform = MatrixRetrieval(inittransform)

    finaltransform = ElastixReg(itk_masked_image,moving_image,fixed_mask,inittransform)
    final_transform = MatrixRetrieval(finaltransform)
    
    Transforms = init_transform + final_transform
    # Transforms = final_transform

    Rotation, Translation = [], []
    for transform in Transforms:
        # print(transform.GetParameters())
        Rotation.append(transform.GetMatrix())
        Translation.append(transform.GetTranslation())
    
    FINAL_TRANSFORM = ComputeFinalMatrix(Rotation,Translation)
    # print(FINAL_TRANSFORM.GetParameters())

    sitk.WriteTransform(FINAL_TRANSFORM,'/home/luciacev/Downloads/TransformITK.tfm')
    # print(matrix_transform)
    # print(transform.GetParameterMap(0))
    # itk.imwrite(im,'/home/luciacev/Desktop/Luc/Projects/Areg_CBCT/im.nii.gz')






def SimpleElastixApprox(fixed_image, moving_image):
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.LogToConsoleOff()
    elastixImageFilter.SetFixedImage(fixed_image)
    elastixImageFilter.SetMovingImage(moving_image)

    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("rigid"))
    elastixImageFilter.SetParameterMap(parameterMapVector)

    # elastixImageFilter.SetParameter("MaximumNumberOfIterations", "5000")
    # elastixImageFilter.SetParameter("NumberOfSpatialSamples", "100000")
    
    tic = time.time()
    elastixImageFilter.Execute()
    # print('Process Time: {} sec'.format(round(time.time()-tic,2)))

    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    return resultImage, transformParameterMap

def SimpleElastixReg(fixed_image, moving_image):
    """Get the transform to register two images using SimpleElastix registration method"""
    
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.LogToConsoleOff()
    elastixImageFilter.SetFixedImage(fixed_image)
    elastixImageFilter.SetMovingImage(moving_image)

    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("rigid"))
    # parameterMapVector.append(sitk.GetDefaultParameterMap("rigid"))
    elastixImageFilter.SetParameterMap(parameterMapVector)
    
    elastixImageFilter.SetParameter("ErodeMask", "true")
    elastixImageFilter.SetParameter("MaximumNumberOfIterations", "10000")
    elastixImageFilter.SetParameter("NumberOfSpatialSamples", "10000")
    elastixImageFilter.SetParameter("NumberOfResolutions", "1")
    
    tic = time.time()
    elastixImageFilter.Execute()
    # print('Process Time: {} sec'.format(round(time.time()-tic,2)))

    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    return resultImage, transformParameterMap

def VoxelBasedRegistration(fixed_image_path,moving_image_path,fixed_seg_path,approx=False):

    # Copy T1 and T2 images to output directory
    # shutil.copyfile(fixed_image_path, os.path.join(outpath,patient+'_T1.nii.gz'))
    # shutil.copyfile(moving_image_path, os.path.join(outpath,patient+'_T2.nii.gz'))
    
    # Read images and segmentations
    fixed_image = sitk.ReadImage(fixed_image_path)
    fixed_seg = sitk.ReadImage(fixed_seg_path)
    fixed_seg.SetOrigin(fixed_image.GetOrigin())
    moving_image = sitk.ReadImage(moving_image_path)
    # moving_seg = sitk.ReadImage(moving_seg_path)

    # Apply mask to images
    fixed_image_masked = applyMask(fixed_image, fixed_seg)
    # moving_image_masked = applyMask(moving_image, moving_seg)

    
    # Register images
    Transforms = []

    if approx:
        # Approximate registration
        # tic = time.time()
        resample_approx, TransformParamMap = SimpleElastixApprox(fixed_image, moving_image)
        Transforms_Approx = MatrixRetrieval(TransformParamMap)
        # print('Registration time: ', round(time.time() - tic,2),'s')
        Transforms = Transforms_Approx
    else:
        resample_approx = moving_image
    
    # Fine tuning
    # tic = time.time()
    resample_t2, TransformParamMap = SimpleElastixReg(fixed_image_masked, resample_approx)
    Transforms_Fine = MatrixRetrieval(TransformParamMap)

    # Combine transforms
    Transforms += Transforms_Fine
    transform = sitk.Transform()
    for t in Transforms:
        transform.AddTransform(t)
    # Resample images and segmentations using the final transform
    # tic = time.time()
    resample_t2 = sitk.Cast(ResampleImage(moving_image, transform),sitk.sitkInt16)
    # resample_t2_seg = ResampleImage(moving_seg, transform)

    
    # Compare segmentations
    # DiceStart,JaccardStart,HausdorfStart = CompareScans(fixed_seg, moving_seg)
    # DiceEnd,JaccardEnd,HausdorfEnd = CompareScans(fixed_seg, resample_t2_seg)
    # print("Delta Dice {} | Delta Jaccard {} | Delta Hausdorf {}".format(round(DiceEnd-DiceStart,5),round(JaccardEnd-JaccardStart,5),round(HausdorfEnd-HausdorfStart,5)))
    # sitk.WriteImage(sitk.Cast(resample_t2,sitk.sitkInt16), os.path.join(outpath,patient+'_ScanReg.nii.gz'))
    # print('Resampling time: ', round(time.time() - tic,2),'s')

    return transform, resample_t2