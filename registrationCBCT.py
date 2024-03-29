import SimpleITK as sitk
import os
import argparse
from tqdm import tqdm

from utils import GetDictPatients, LoadOnlyLandmarks, applyTransformLandmarks, WriteJson, VoxelBasedRegistration, translate

def main(args):
    t1_folder,t2_folder, output_dir = args.t1_folder, args.t2_folder, args.output_dir

    patients = GetDictPatients(t1_folder,t2_folder,args.reg_type,args.todo)
    
    for i,(patient,data) in tqdm(enumerate(patients.items()),total=len(patients)):
        if args.print:
            print("="*70)
            print("T1 scan path: {} \nT2 scan path: {}".format(os.path.basename(data['scanT1']),os.path.basename(data['scanT2'])))
            print("T1 seg path: {}".format(os.path.basename(data['segT1'])))
            if args.reg_lm:
                print("T2 landmarks path: {}".format(os.path.basename(data['lmT2'])))
        # print("Working on patient {}".format(patient))
        else:
            try:
                transform, resample_t2 = VoxelBasedRegistration(data['scanT1'],data['scanT2'],data['segT1'],args.approx)
                
                outpath = os.path.join(output_dir,translate(args.reg_type),patient+'_OutReg')
                if not os.path.exists(outpath):
                    os.makedirs(outpath)
                sitk.WriteTransform(transform, os.path.join(outpath,patient+'_matrix.tfm'))
                sitk.WriteImage(resample_t2, os.path.join(outpath,patient+'_ScanReg.nii.gz'))
                if args.reg_lm:   
                    transformedLandmarks = applyTransformLandmarks(LoadOnlyLandmarks(data['lmT2']), transform.GetInverse())
                    WriteJson(transformedLandmarks, os.path.join(outpath,patient+'_lm_Reg.mrk.json'))
            
            except KeyError:
                # print("Patient {} does not have both T1 and T2 scans".format(patient))
                continue

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Voxel Based Registration')
    
    parser.add_argument('--t1_folder', type=str, help='Path to folder containing input T1 scans',default='/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/JJ/Approach3/T1Or/')
    parser.add_argument('--t2_folder', type=str, help='Path to folder containing input T2 scans',default='/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/JJ/Approach3/T2/')
    parser.add_argument('--output_dir', type=str, help='Path to folder containing output register T2 scans',default='/home/luciacev/Desktop/Luc/DATA/AReg_CBCT/NOT WORKING/')
    parser.add_argument("--reg_type", type=str, help="Type of registration to perform", default="MAND", choices=['CB','MAND','MAX'])
    parser.add_argument("--print", type=bool, help="Print info", default=False)
    parser.add_argument("--todo", type=str, help="What scan to do", default='')
    parser.add_argument("--reg_lm", type=bool, help="Whether or not we apply the matrix to the landmark file", default=False)
    parser.add_argument("--approx", type=bool, help="Whether or not we use the approximate method", default=False)
    
    args = parser.parse_args()

    main(args)