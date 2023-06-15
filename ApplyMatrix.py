import SimpleITK as sitk
from utils import WriteJson, LoadOnlyLandmarks, applyTransformLandmarks

transform_path = '/home/luciacev/Desktop/Luc/DATA/SOPHIE/14/TB_0014_matrix.tfm'
transform = sitk.ReadTransform(transform_path)
json_path = '/home/luciacev/Desktop/Luc/DATA/SOPHIE/14/TB_0014_T2_lm_Or.mrk.json'
transformedLandmarks = applyTransformLandmarks(LoadOnlyLandmarks(json_path), transform.GetInverse())
WriteJson(transformedLandmarks, '/home/luciacev/Desktop/Luc/DATA/SOPHIE/14/TB_0014_T2_lm_Reg.mrk.json')