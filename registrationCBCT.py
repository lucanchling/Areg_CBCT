import SimpleITK as sitk
import os
import argparse
from tqdm import tqdm

from utils import *

def main(args):
    t1_folder,t2_folder, output_dir = args.t1_folder, args.t2_folder, args.output_dir

    t1_patients = GetPatients(t1_folder, time_point='T1',segmentationType=args.reg_type)
    t2_patients = GetPatients(t2_folder, time_point='T2',segmentationType=args.reg_type)

    # merge dictionaries for each patient
    patients = MergeDicts(t1_patients,t2_patients)

    # for patient,data in patients.items():
    for i,(patient,data) in tqdm(enumerate(patients.items()),total=len(patients)):
        # print("="*70)
        # print("Working on patient {}".format(patient))
        try:
            # print("T1 scan path: {} \nT2 scan path: {}".format(os.path.basename(data['scanT1']),os.path.basename(data['scanT2'])))
            # print("T1 seg path: {} \nT2 seg path: {}".format(os.path.basename(data['segT1']),os.path.basename(data['segT2'])))
            # print("T2 landmarks path: {}".format(os.path.basename(data['lmT2'])))
            outpath = os.path.join(output_dir,translate(args.reg_type),patient+'_OutReg')
            transform, resample_t2, resample_t2_seg = VoxelBasedRegistration(data['scanT1'],data['scanT2'],data['segT1'],data['segT2'],outpath,patient)
            transformedLandmarks = applyTransformLandmarks(LoadOnlyLandmarks(data['lmT2']), transform.GetInverse())
            sitk.WriteImage(resample_t2, os.path.join(outpath,patient+'_ScanReg.nii.gz'))
            sitk.WriteImage(resample_t2_seg, os.path.join(outpath,patient+'_'+args.reg_type+'MASK_SegReg.nii.gz'))
            WriteJson(transformedLandmarks, os.path.join(outpath,patient+'_lm_Reg.mrk.json'))
        except KeyError:
            print("Patient {} does not have both T1 and T2 scans".format(patient))
            continue

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Voxel Based Registration')
    parser.add_argument('--t1_folder', type=str, help='Path to folder containing input T1 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/TwinBlock/3_TB1Or/')
    parser.add_argument('--t2_folder', type=str, help='Path to folder containing input T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/TwinBlock/3_TB1.5Or/')
    parser.add_argument('--output_dir', type=str, help='Path to folder containing output register T2 scans',default='/home/lucia/Desktop/Luc/DATA/AReg/SOPHIE/TEST/OUT/')
    parser.add_argument("--reg_type", type=str, help="Type of registration to perform", default='CB', choices=['CB','MAND','MAX'])

    args = parser.parse_args()

    main(args)