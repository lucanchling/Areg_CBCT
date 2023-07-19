import SimpleITK as sitk
import numpy as np

def MatrixRetrieval(transform_file):
    """Retrieve the matrix from the transform file"""

    Rotation = []
    Translation = []

    # Read the transform file as a txt file
    with open(transform_file, 'r') as file:
        lines = file.readlines()[1:]
        # Clean the \n at the end of each line
        lines = [line.strip() for line in lines]
        for i,line in enumerate(lines):
            if line.startswith('Transform:'):
                transform_type = line.split(' ')[1]
                if transform_type == 'Euler3DTransform_double_3_3':
                    parameters = [float(i) for i in lines[i+1].split('Parameters: ')[1].split(' ')]
                    A = parameters[0:3]
                    B = parameters[3:6]
                    # Convert to a sitk transform
                    transform = sitk.Euler3DTransform()
                    transform.SetRotation(angleX=A[0], angleY=A[1], angleZ=A[2])
                    transform.SetTranslation(B)
                    Rotation.append(transform.GetMatrix())
                    Translation.append(transform.GetTranslation())
                    # print('Matrix:',matrix)
                    # print('Translation:',translation)
    
    return Rotation, Translation

def ComputeFinalMatrix(Rotation,Translation):
    """Compute the final matrix from the list of matrices and translations"""

    # Compute the final rotation matrix
    final_rotation = np.reshape(np.asarray(Rotation[0]),(3,3))
    for i in range(1,len(Rotation)):
        final_rotation = final_rotation @ np.reshape(np.asarray(Rotation[i]),(3,3))
    
    # Compute the final translation matrix
    final_translation = np.reshape(np.asarray(Translation[0]),(1,3))
    for i in range(1,len(Translation)):
        final_translation = final_translation + np.reshape(np.asarray(Translation[i]),(1,3))

    # Create the final transform
    final_transform = sitk.Euler3DTransform()
    final_transform.SetMatrix(final_rotation.flatten().tolist())
    final_transform.SetTranslation(final_translation[0].tolist())

    return final_transform

def GetAngleFromRotationMatrix(Rotation):
    """Compute the angles from the rotation matrix"""

    # Compute the angles
    angle_x = np.arctan2(Rotation[2,1],Rotation[2,2])
    angle_y = np.arctan2(-Rotation[2,0],np.sqrt(Rotation[2,1]**2 + Rotation[2,2]**2))
    angle_z = np.arctan2(Rotation[1,0],Rotation[0,0])

    return np.array([angle_x,angle_y,angle_z])

def main(args):

    Rot, Trans = MatrixRetrieval(args.transform)

    fin_Rot, fin_Trans = ComputeFinalMatrix(Rot,Trans)

    # print('Final Rotation Matrix:\n',fin_Rot)
    # print('Final Translation Matrix:\n',fin_Trans)

    angles = GetAngleFromRotationMatrix(fin_Rot)
    # print('Angles:',angles)

    # Compute the final transform
    final_transform = sitk.Euler3DTransform()
    final_transform.SetRotation(angles[0],angles[1],angles[2])
    final_transform.SetTranslation(fin_Trans[0])

    # Save the final transform
    sitk.WriteTransform(final_transform,args.transform)
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--transform', type=str, help='Path to the transform file')
    args = parser.parse_args()

    main(args)




    