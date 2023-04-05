import SimpleITK as sitk
import os
import numpy as np
from utils import GetDictPatients
import argparse
from tqdm import tqdm

def CloseCBCTSeg(input_img, closing_radius = 15):
    """
    Close the holes in the CBCT

    Parameters
    ----------
    filePath
     path of the image file 
    radius
     radius of the closing to apply to the seg
    outpath
     path to save the new image
    """
    img = sitk.GetArrayFromImage(input_img)

    img = np.where(img > 0, 1,img)
    output = sitk.GetImageFromArray(img)
    output.SetSpacing(input_img.GetSpacing())
    output.SetDirection(input_img.GetDirection())
    output.SetOrigin(input_img.GetOrigin())

    output = sitk.BinaryDilate(output, [closing_radius] * output.GetDimension())
    output = sitk.BinaryFillhole(output)
    output = sitk.BinaryErode(output, [closing_radius] * output.GetDimension())

    return output

def CreateMask(args,patients):
    for patient, data in tqdm(patients.items()):
        
        # t1_scan = sitk.ReadImage(data['scanT1'])
        full_seg_t1 = sitk.ReadImage(data['segT1'])
        target_mask_t2 = sitk.ReadImage(data['segT2'])
        Transform = sitk.ReadTransform(data['mat'])

        # Resample the seg to the T1
        resample_seg = sitk.Resample(target_mask_t2,full_seg_t1, Transform, sitk.sitkLinear, 0.0, full_seg_t1.GetPixelID())

        # Use the T2 mask as a filter to create the T1 mask
        crop = sitk.Mask(full_seg_t1, resample_seg)

        # crop.SetSpacing(t1_scan.GetSpacing())
        # crop.SetDirection(t1_scan.GetDirection())
        # crop.SetOrigin(t1_scan.GetOrigin())
        OutPath = os.path.join(os.path.dirname(data['segT1']),'{}_{}MASK_Seg.nii.gz'.format(patient,args.reg_type))
        sitk.WriteImage(crop, OutPath)
        

def main(args):
    data_dir = args.data_dir
    t1_folder,t2_folder,matrix_folder = os.path.join(data_dir,args.reg_type,'T1'), os.path.join(data_dir,args.reg_type,'T2'), os.path.join(data_dir,args.reg_type,'Matrix')

    patients = GetDictPatients(folder_t1_path=t1_folder, folder_t2_path=t2_folder,matrix_folder=matrix_folder)

    CreateMask(args,patients)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create mask for the CBCT')
    parser.add_argument('--data_dir', type=str, help='Folder with the matrix',default='/home/luciacev/Desktop/Luc_Anchling/DATA/AReg_CBCT/JJ/Approach2')
    parser.add_argument('--reg_type', type=str, help='Type of registration',choices=['MAND','CB','MAX'],default='MAND')
    args = parser.parse_args()
    main(args)
