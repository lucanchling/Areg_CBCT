import slicer
import sys


def aff():
    print()
    print("========================================")
    print()


fixedVolume = slicer.util.loadVolume(
    "/home/lucia/Desktop/Luc/DATA/AReg/TestFiles/T1/T1_scan.gipl.gz"
)
movingVolume = slicer.util.loadVolume(
    "/home/lucia/Desktop/Luc/DATA/AReg/TestFiles/T2/T2_Mand_scan.gipl.gz"
)
fixedMaskVolume = slicer.util.loadVolume(
    "/home/lucia/Desktop/Luc/DATA/AReg/TestFiles/T1/T1_Pred_Mand.gipl.gz"
)
movingMaskVolume = slicer.util.loadVolume(
    "/home/lucia/Desktop/Luc/DATA/AReg/TestFiles/T2/T2_Pred_Mand.gipl.gz"
)

outputTransform = slicer.vtkMRMLLinearTransformNode()
slicer.mrmlScene.AddNode(outputTransform)

outputVolume = slicer.vtkMRMLScalarVolumeNode()
slicer.mrmlScene.AddNode(outputVolume)

outputSegmentation = slicer.vtkMRMLScalarVolumeNode()
slicer.mrmlScene.AddNode(outputSegmentation)


cli = {
    "growing": {
        "module": slicer.modules.growing,
        "parameters": {
            "transformPath": outputTransform,
            "movingVolume": movingVolume,
            "fixedVolume": fixedVolume,
            "fixedMaskVolume": fixedMaskVolume,
            "movingMaskVolume": movingMaskVolume,
            "segmentation": movingMaskVolume,
            "segmentationOut": outputSegmentation,
        },
    },
    "brainsfit": {
        "module": slicer.modules.brainsfit,
        "parameters": {
            "fixedVolume": fixedVolume,
            "movingVolume": movingVolume,
            "fixedMaskVolume": fixedMaskVolume,
            "movingMaskVolume": movingMaskVolume,
            "useRigid": True,
            "linearTransform": outputTransform,
            "outputVolume": outputVolume,
        },
    },
}
module = "brainsfit"
node = slicer.cli.run(cli[module]["module"], parameters=cli[module]["parameters"])


aff()
quit()
