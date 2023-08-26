import numpy as np
from Femur_MA import femur_MA
from OpenSimMarkers import OpenSimMarkers
from coordinatesCorrection import coordinatesCorrection
from femurShaft_ns import femurShaft_ns
import matplotlib.pyplot as plt
from struct2xml import struct2xml
from coordinatesOpenSim import coordinatesOpenSim
import os
from scipy.spatial import Delaunay

def femur_ns(dataModel, markerset, answerLeg, rightbone, FA_angle, NS_angle, answerNameModelFemur,
             answerNameMarkerFemur, dataFemur, place):
			 
	# Bone vertices
    # Convert the vertices into numpy array and find the polygons


    femur = np.array(dataFemur["PolyData"]["Piece"]["Points"]["DataArray"].split()).astype(float)
    femur_2 = [femur[i:i+3] for i in range(0, len(femur), 3)]
    polyText = dataFemur["PolyData"]["Piece"]["Polys"]["DataArray"][0]
    
    # print(polyText[1])
    # The muscle attachments for the femur are put in one matrix.
    femurMuscle, femurPlace1, femurNR, femurMuscleType = femur_MA(dataModel, answerLeg, rightbone)
    # print(femurMuscle)
    # Find the markers attached to the femur in OpenSim
    _, _, _, _, markerFemur, markerFemurNR = OpenSimMarkers(markerset, answerLeg, rightbone)
    
    
    
    # % The vertices for the bone and muscle attachements
    # % are rotated to fit the coordinate system in MATLAB
    femur_start = coordinatesCorrection(femur_2)
    femurMuscle_start = coordinatesCorrection(femurMuscle)
    # print('muscle start')
    # print(femurMuscle_start.T)
    markerFemur_start = coordinatesCorrection(markerFemur)

    wrapCnt = 0
    wrapLocations = []
    wrapRotations = []
    wrapIndizes = []
        


    for i in range(1,len(dataModel["Model"]["BodySet"]["objects"]["Body"])):
        if 'femur_' + answerLeg.lower() in dataModel["Model"]["BodySet"]["objects"]["Body"][i]["Joint"]["CustomJoint"]["parent_body"]:
            #print("yes 1")
            if  "WrapObjectSet" in dataModel["Model"]["BodySet"]["objects"]["Body"][i]:
                #print("yes 2")
                if "objects" in dataModel["Model"]["BodySet"]["objects"]["Body"][i]["WrapObjectSet"]:
                    #print("yes 3")
                    wrapObjectTypes = dataModel["Model"]["BodySet"]["objects"]["Body"][i]["WrapObjectSet"]["objects"]["WrapCylinder"]
                    for j in range(len(wrapObjectTypes)):
                        wrapRotations.append(list(map(float, dataModel["Model"]["BodySet"]["objects"]["Body"][i]["WrapObjectSet"]["objects"]["WrapCylinder"][j]["xyz_body_rotation"].split())))
                        wrapLocations.append(list(map(float, dataModel["Model"]["BodySet"]["objects"]["Body"][i]["WrapObjectSet"]["objects"]["WrapCylinder"][j]["xyz_body_rotation"].split())))
                        for k in range(len(dataModel["Model"]["BodySet"]["objects"]["Body"][i]["WrapObjectSet"]["objects"]["WrapCylinder"][j])):
                            wrapCnt += 1
                            wrapIndizes.append([i, j, k])
            break
    # print(wrapLocations)
    # print(wrapRotations)
    # print(wrapIndizes)

    wrapRotations_New = wrapRotations.copy()
    wrapLocations_New = wrapLocations.copy()
    wrapLocations = coordinatesCorrection(wrapLocations)

    innerBox, middleBox, innerBoxMA, innerBoxMarker, middleBoxMA, middleBoxMarker, femur_NewAxis, H_transfer, \
    angleZX, angleZY, angleXY, centroidValueLGtroch, femurMA_NewAxis, femurMarker_NewAxis, \
    Condyl_NewAxis, Shaft_proximal, Shaft_distal, CondylMA_NewAxis, ShaftMA_distal, CondylMarker_NewAxis, ShaftMarker_distal, \
    wrapLocations_NewAxis = femurShaft_ns(dataModel, femur_start, answerLeg, rightbone, femurMuscle_start, markerFemur_start, wrapLocations)
    #print(Condyl_NewAxis[0][0])   
    
    ## import and convert the polies for figures
    femurSize = femur_start.shape
    polysplit = polyText.split('\n')
    
    #print(polysplit)
    poly3 = []
    for i in range(len(polysplit)):
        poly3.append(np.array(polysplit[i].split()).astype(float))
    polys = []
    for i in range(1,len(poly3) - 1):
        polys.append(poly3[i] + 1)
    polys = np.array(polys)
    # print(polys)
    # polys = np.array(polys)
    
    ## Rotation; step 1
    
    Ry_NS = np.array([[np.cos(NS_angle), 0, np.sin(NS_angle)],
                  [0, 1, 0],
                  [-np.sin(NS_angle), 0, np.cos(NS_angle)]])  # neck shaft angle
    Rz_FA = np.array([[np.cos(FA_angle), -np.sin(FA_angle), 0],
                      [np.sin(FA_angle), np.cos(FA_angle), 0],
                      [0, 0, 1]])  # femoral anteversion
    RotMatrix = np.dot(Rz_FA, Ry_NS)

    femur_rot1_all = []
    innerBox_rot1 = []
    polys_innerNumber = []
    for i in range(len(femur_NewAxis)):
        if np.any(np.all(femur_NewAxis[i, :] == innerBox, axis=1)):
            item_innerBox = np.dot(RotMatrix, femur_NewAxis[i, :].reshape(-1, 1)).flatten()
            femur_rot1_all.append(item_innerBox)
            innerBox_rot1.append(item_innerBox)
            polys_innerNumber.append(i + 1)
        else:
            femur_rot1_all.append(femur_NewAxis[i, :])
            
    
    femur_rot1_all_2 = femur_rot1_all  
    #print(femur_rot1_all_2)
    femur_rot1_all = np.array(femur_rot1_all)
    #print(femur_rot1_all)
    
# =============================================================================
#     const chunkSize = 1;
#     for (let i = 0; i < array.length; i += chunkSize) {
#     const chunk = array.slice(i, i + chunkSize);
#     // do whatever
#     }
# =============================================================================
    
    
    # print(femur_rot1_all)
    # for i in femur_NewAxis:
    #     print(i)
    # for i in femur_rot1_all: 
    #     print(i)

    # zeroMatrix = np.zeros_like(polys)
    # for ii in range(len(polys_innerNumber)):
    #     zeroMatrix += (polys == polys_innerNumber[ii])

    # polys_inner = []
    # for k in range(zeroMatrix.shape[0]):
    #     if np.sum(zeroMatrix[k, :]) == 3:
    #         polys_inner.append(polys[k, :])
    # polys_inner = np.array(polys_inner)

    femurMA_rot1_all = []
    for i in range(len(femurMA_NewAxis)):
        if np.any(np.all(femurMA_NewAxis[i, :] == innerBoxMA, axis=1)):
            item_innerBoxMA = np.dot(RotMatrix, femurMA_NewAxis[i, :].reshape(-1, 1)).flatten()
            femurMA_rot1_all.append(item_innerBoxMA)
        else:
            femurMA_rot1_all.append(femurMA_NewAxis[i, :])
    femurMA_rot1_all = np.array(femurMA_rot1_all)

    femurMarker_rot1_all = []
    for i in range(len(femurMarker_NewAxis)):
        if np.any(np.all(femurMarker_NewAxis[i, :] == innerBoxMarker, axis=1)):
            item_innerBoxMarker = np.dot(RotMatrix, femurMarker_NewAxis[i, :].reshape(-1, 1)).flatten()
            femurMarker_rot1_all.append(item_innerBoxMarker)
        else:
            femurMarker_rot1_all.append(femurMarker_NewAxis[i, :])
    femurMarker_rot1_all = np.array(femurMarker_rot1_all)

    innerBox_H_rot = np.dot(RotMatrix, H_transfer.reshape(-1, 1)).flatten()
    
# =============================================================================
    def set_axes_equal(ax):
    # """
    # Make axes of 3D plot have equal scale so that spheres appear as spheres,
    # cubes as cubes, etc.

    # Input
    #   ax: a matplotlib axis, e.g., as output from plt.gca().
    # """

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()
    
        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)
    
        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])
    
        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
# =============================================================================
    
    fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111, projection='3d')    #original
    ax = fig.add_subplot(111, projection='3d')
    # ax = plt.axes(projection='3d')
    # triangles=Delaunay(femur_rot1_all[:,:2]).simplices
    ax.scatter(femur_rot1_all[:, 0], femur_rot1_all[:, 1], femur_rot1_all[:, 2], edgecolor='black', linestyle=':')
    ax.scatter(femur_NewAxis[:, 0], femur_NewAxis[:, 1], femur_NewAxis[:, 2], edgecolor='black')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_box_aspect([1.0, 1.0, 1.0])
    set_axes_equal(ax)
    
    maxZ = np.max(middleBox[:, 2])
    minZ = np.min(middleBox[:, 2])

    # #     for k in range(wrapCnt):
    # #         print(k)


    for u in range(wrapCnt):
        if u < 5:
            currLoc = wrapLocations_NewAxis[u, :]
            if currLoc[2] > minZ and currLoc[2] < maxZ:
                i = wrapIndizes[u, 0]
                j = wrapIndizes[u, 1]
                k = wrapIndizes[u, 2]
                name = dataModel["Model"]["BodySet"]["objects"]["Body"][0, i]["WrapObjectSet"]["objects"].wrapObjectTypes[j][k].Attributes.name
                print('WrapObject "' + name + '" is in the area where the femur itself is rotated! This is not supported yet, adjust the location manually!')
    
    ## Rotation: step 2
    
    high_middleBox = np.max(middleBox[:, 2])
    low_middleBox = np.min(middleBox[:, 2])
    low_value = np.zeros(3)
    high_value = np.zeros(3)
    
    for i in range(middleBox.shape[0]):
        if middleBox[i, 1] <= low_middleBox:
            low_value = middleBox[i, :]
        else:
            high_value = middleBox[i, :]
    
    topmiddlebox = np.array([0, 0, high_middleBox])
    bottommiddlebox = np.array([0, 0, low_middleBox])
    middlevector = bottommiddlebox - topmiddlebox
    distanceMiddleBox = np.linalg.norm(middlevector)
    polys_middle_Number = []
    femur_rot2_all = np.zeros_like(femur_rot1_all)
    middleBox_rot2 = []
    middleBox_before_rot = []
    
    for i in range(femur_rot1_all.shape[0]):
        if np.any(np.all(femur_rot1_all[i] == middleBox, axis=1)):
            middleBox_before_rot.append(femur_rot1_all[i])
            a_proj = femur_rot1_all[i][2] - topmiddlebox
            projv = np.dot(a_proj, middlevector) / distanceMiddleBox
            scalingRotVectZ = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox) * FA_angle
            RotMatrix_rot2 = np.array([[np.cos(scalingRotVectZ), -np.sin(scalingRotVectZ), 0],
                                       [np.sin(scalingRotVectZ), np.cos(scalingRotVectZ), 0],
                                       [0, 0, 1]])
            scalingRotVectY = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox)
            grad_rot2_middleBox = np.dot(RotMatrix_rot2, femur_rot1_all[i, :])
            femur_rot2_all[i, :] = grad_rot2_middleBox
            middleBox_rot2.append(grad_rot2_middleBox)
            polys_middle_Number.append(i)
        else:
            femur_rot2_all[i, :] = femur_rot1_all[i, :]
    
    zeroMatrix = np.zeros_like(polys)
    for ii in polys_middle_Number:
        zeroMatrix += (polys == ii)
    
    polys_middle = []
    for k in range(zeroMatrix.shape[0]):
        if np.sum(zeroMatrix[k, :]) == 3:
            polys_middle.append(polys[k, :])
    
    femurMA_rot2_all = np.zeros_like(femurMA_rot1_all)
    middleBoxMA_rot2 = []
    middleBoxMA_before_rot = []
    
    for i in range(femurMA_rot1_all.shape[0]):
        if np.any(np.all(femurMA_rot1_all[i, :] == middleBoxMA, axis=1)):
            middleBoxMA_before_rot.append(femurMA_rot1_all[i])
            a_proj = femurMA_rot1_all[i, :3] - topmiddlebox
            projv = np.dot(a_proj, middlevector) / distanceMiddleBox
            scalingRotVectZ = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox) * FA_angle
            RotMatrix_rot2 = np.array([[np.cos(scalingRotVectZ), -np.sin(scalingRotVectZ), 0],
                                       [np.sin(scalingRotVectZ), np.cos(scalingRotVectZ), 0],
                                       [0, 0, 1]])
            scalingRotVectY = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox)
            grad_rot2_middleBoxMA = np.dot(RotMatrix_rot2, femurMA_rot1_all[i, :])
            femurMA_rot2_all[i, :] = grad_rot2_middleBoxMA
            middleBoxMA_rot2.append(grad_rot2_middleBoxMA)
        else:
            femurMA_rot2_all[i, :] = femurMA_rot1_all[i, :]
    
    femurMarker_rot2_all = np.zeros_like(femurMarker_rot1_all)
    middleBoxMarker_rot2 = []
    middleBoxMarker_before_rot = []
    print(middleBoxMarker)
    for i in range(len(femurMarker_rot1_all)):
        if np.any(np.all(femurMarker_rot1_all[i] == middleBoxMarker, axis=1)):
            middleBoxMarker_before_rot = np.vstack((middleBoxMarker_before_rot, femurMarker_rot1_all[i, :]))
            a_proj = femurMarker_rot1_all[i, :3] - topmiddlebox
            projv = np.dot(a_proj, middlevector) / distanceMiddleBox
            scalingRotVectZ = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox) * FA_angle
            RotMatrix_rot2 = np.array([[np.cos(scalingRotVectZ), -np.sin(scalingRotVectZ), 0],
                                       [np.sin(scalingRotVectZ), np.cos(scalingRotVectZ), 0],
                                       [0, 0, 1]])
            scalingRotVectY = (np.abs(projv - distanceMiddleBox) / distanceMiddleBox)
            grad_rot2_middleBoxMarker = np.dot(RotMatrix_rot2, femurMarker_rot1_all[i, :])
            femurMarker_rot2_all[i, :] = grad_rot2_middleBoxMarker
            middleBoxMarker_rot2 = np.vstack((middleBoxMarker_rot2, grad_rot2_middleBoxMarker))
        else:
            femurMarker_rot2_all[i, :] = femurMarker_rot1_all[i, :]
    
    fig = plt.figure(figsize=(10, 20))
    ax1 = fig.add_subplot(211, projection='3d')
    ax1.plot_trisurf(polys, femur_rot2_all[:, 0], femur_rot2_all[:, 1], femur_rot2_all[:, 2], color='none', edgecolor='black', linestyle=':')
    ax1.plot_trisurf(polys, femur_NewAxis[:, 0], femur_NewAxis[:, 1], femur_NewAxis[:, 2], color='none', edgecolor='black')
    ax1.set_box_aspect([1, 1, 1])
    ax1.view_init(30, -10)
    ax1.grid(True)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_zticklabels([])
    
    ax2 = fig.add_axes([0.65, 0.25, 0.35, 0.55], projection='3d')
    ax2.plot_trisurf(polys_middle, femur_rot2_all[:, 0], femur_rot2_all[:, 1], femur_rot2_all[:, 2], color='none', edgecolor='black', linestyle=':')
    ax2.plot_trisurf(polys_middle, femur_NewAxis[:, 0], femur_NewAxis[:, 1], femur_NewAxis[:, 2], color='none', edgecolor='black')
    ax2.plot_trisurf(polys_inner, femur_rot2_all[:, 0], femur_rot2_all[:, 1], femur_rot2_all[:, 2], color='none', edgecolor='black', linestyle=':')
    ax2.plot_trisurf(polys_inner, femur_NewAxis[:, 0], femur_NewAxis[:, 1], femur_NewAxis[:, 2], color='none', edgecolor='black')
    ax2.set_box_aspect([1, 1, 1])
    ax2.view_init(30, -10)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_zticklabels([])
    
    maxZ = np.max(middleBox[:, 2])
    minZ = np.min(middleBox[:, 2])
    for u in range(wrapCnt):
        currLoc = wrapLocations_NewAxis[u, :]
        if minZ < currLoc[2] < maxZ:
            i = wrapIndizes[u, 0]
            j = wrapIndizes[u, 1]
            k = wrapIndizes[u, 2]
            name = dataModel.OpenSimDocument.Model.BodySet.objects.Body[0, i].WrapObjectSet.objects[wrapObjectTypes[j]][0, k].Attributes.name
            print(f'WrapObject "{name}" is in the area where the femur itself is rotated! This is not supported yet, adjust the location manually!')
            
            
    # Rotation; Step 3
    # %MOVE THE FEMORAL HEAD BACK TO FIT CONDYLAR WITHOUT MOVING THE CONDYLAR
    # %start by transfering the bone back to the H_ZX
    
    transfer_step3 = H_transfer - innerBox_H_rot

    # find the distance between the highest value in the distal shaft to the highest value in the condylar
    maxCondyl = np.max(Condyl_NewAxis[:][2])
    condyl_top = np.zeros(3)
    for i in range(len(Condyl_NewAxis)):
        if Condyl_NewAxis[i][2] >= maxCondyl:
            condyl_top = Condyl_NewAxis[i]
            
    minShaft_prox = np.min(Shaft_proximal[:][2])
    shaft_prox_min = np.zeros(3)
    for i in range(len(Shaft_proximal)):
        if Shaft_proximal[i][2] <= minShaft_prox:
            shaft_prox_min = Shaft_proximal[i]

    distance_shaft_distal = np.linalg.norm(condyl_top - shaft_prox_min)

    femur_rot3_all = []  # The outer box is moved back to restore the position of the femoral head to the acetabulum
    femur_rot3_all_deform = []  # the outer box has been restored and the distal part of the femoral shaft is gradually deformed to fit the condylar (which does not move)
    
    
    for i in range(femur_rot2_all.shape[0]):
        if np.any(np.all(femur_rot2_all[i] == Condyl_NewAxis, axis=1)):  # the condylar do not move but the rest is translated to restore the femoral head
            femur_rot3_all.append(femur_rot2_all[i])
        else:
            item_rot3 = femur_rot2_all[i] + transfer_step3  # The outer box is transferred back to the position of the femoral head
            femur_rot3_all.append(item_rot3)

        if np.any(np.all(femur_rot2_all[i] == Shaft_distal, axis=1)):
            scaler = (np.linalg.norm(femur_rot3_all[i] - condyl_top) - distance_shaft_distal) / distance_shaft_distal * transfer_step3
            item_rot3_distal = femur_rot3_all[i] + scaler
            femur_rot3_all_deform.append(item_rot3_distal)
        else:
            # print(femur_rot3_all)  
            # print(femur_rot3_all_deform)
            femur_rot3_all_deform.append(femur_rot3_all[i])
            
    
    femurMA_rot3_all = np.zeros_like(femurMA_rot2_all)
    femurMA_rot3_all_deform = np.zeros_like(femurMA_rot2_all)
    test = []

    for i in range(femurMA_rot2_all.shape[0]):
        if np.any(np.all(femurMA_rot2_all[i, :] == CondylMA_NewAxis, axis=1)):  # the condylar do not move
            femurMA_rot3_all[i, :] = femurMA_rot2_all[i, :]
            test.append(femurMA_rot2_all[i, :])
        else:
            item_rot3 = femurMA_rot2_all[i, :] + transfer_step3  # The inner and middle box are transferred back to the position of the femoral head
            femurMA_rot3_all[i, :] = item_rot3

        if np.any(np.all(femurMA_rot2_all[i, :] == ShaftMA_distal, axis=1)):
            scaler = (np.linalg.norm(femurMA_rot3_all[i, :] - condyl_top) - distance_shaft_distal) / distance_shaft_distal * transfer_step3
            item_rot3_distal = femurMA_rot3_all[i, :] + scaler
            femurMA_rot3_all_deform[i, :] = item_rot3_distal
        else:
            femurMA_rot3_all_deform[i, :] = femurMA_rot3_all[i, :]
            
            
    femurMarker_rot3_all = np.zeros_like(femurMarker_rot2_all)
    femurMarker_rot3_all_deform = np.zeros_like(femurMarker_rot2_all)

    for i in range(femurMarker_rot2_all.shape[0]):
        if np.any(np.all(femurMarker_rot2_all[i, :] == CondylMarker_NewAxis, axis=1)):  # the condylar do not move
            femurMarker_rot3_all[i, :] = femurMarker_rot2_all[i, :]
        else:
            item_rot3 = femurMarker_rot2_all[i, :] + transfer_step3  # The inner and middle box are transferred back to the position of the femoral head
            femurMarker_rot3_all[i, :] = item_rot3

        if np.any(np.all(femurMarker_rot2_all[i, :] == ShaftMarker_distal, axis=1)):
            scaler = (np.linalg.norm(femurMarker_rot3_all[i, :] - condyl_top) - distance_shaft_distal) / distance_shaft_distal * transfer_step3
            item_rot3_distal = femurMarker_rot3_all[i, :] + scaler
            femurMarker_rot3_all_deform[i, :] = item_rot3_distal
        else:
            femurMarker_rot3_all_deform[i, :] = femurMarker_rot3_all[i, :]
            
            
    maxShaft_prox = np.max(Shaft_proximal[:, 2])

    for i in range(wrapCnt):
        currLoc = wrapLocations_NewAxis[i, :]
        wrapRotations_New[i, :] = wrapRotations_New[i, :]

        if currLoc[2] < maxCondyl:  # it is at the condyles, don't do any rotation
            wrapLocations_New[i, :] = wrapLocations_New[i, :]
        else:
            wrapLocations_New[i, :] = wrapLocations_New[i, :] + coordinatesOpenSim(transfer_step3)

        if currLoc[2] < maxShaft_prox and currLoc[2] > minShaft_prox:
            scaler = (np.abs(np.linalg.norm(wrapLocations_New[i, :] - condyl_top)) - distance_shaft_distal) / distance_shaft_distal * transfer_step3
            wrapLocations_New[i, :] = wrapLocations_New[i, :] + coordinatesOpenSim(scaler)
            
            
    for u in range(wrapCnt):
        i = wrapIndizes[u, 0]
        j = wrapIndizes[u, 1]
        k = wrapIndizes[u, 2]
        dataModel["Model"]["BodySet"]["objects"]["Body"][0, i]["WrapObjectSet"]["objects"].wrapObjectTypes[j][0, k].xyz_body_rotation.Text = str(wrapRotations_New[u, :])
        dataModel["Model"]["BodySet"]["objects"]["Body"][0, i]["WrapObjectSet"]["objects"].wrapObjectTypes[j][0, k].translation.Text = str(wrapLocations_New[u, :])
        
        
    if femurMA_NewAxis.size != 0:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(femurMA_NewAxis[:, 0], femurMA_NewAxis[:, 1], femurMA_NewAxis[:, 2], c='b', s=20)
        ax.scatter(femurMA_rot3_all_deform[:, 0], femurMA_rot3_all_deform[:, 1], femurMA_rot3_all_deform[:, 2], c='black', s=20)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()


    angleZX_back = -angleZX
    angleZY_back = -angleZY

    Rx_FA = np.array([[1, 0, 0],
                      [0, np.cos(angleZY_back), -np.sin(angleZY_back)],
                      [0, np.sin(angleZY_back), np.cos(angleZY_back)]])

    if answerLeg == rightbone:
        Ry_FA = np.array([[np.cos(angleZX_back), 0, np.sin(angleZX_back)],
                          [0, 1, 0],
                          [-np.sin(angleZX_back), 0, np.cos(angleZX_back)]])
    else:
        Ry_FA = np.array([[np.cos(-angleZX_back), 0, np.sin(-angleZX_back)],
                          [0, 1, 0],
                          [-np.sin(-angleZX_back), 0, np.cos(-angleZX_back)]])

    angleXY_back = -angleXY

    if answerLeg == rightbone:
        Rz_FA = np.array([[np.cos(angleXY_back), -np.sin(angleXY_back), 0],
                          [np.sin(angleXY_back), np.cos(angleXY_back), 0],
                          [0, 0, 1]])
    else:
        Rz_FA = np.array([[np.cos(-angleXY_back), -np.sin(-angleXY_back), 0],
                          [np.sin(-angleXY_back), np.cos(-angleXY_back), 0],
                          [0, 0, 1]])

    R_backOpenSim = Rx_FA @ Ry_FA @ Rz_FA

    femur_rot_back = np.zeros_like(femur_rot3_all_deform)
    for i in range(femur_rot3_all_deform.shape[0]):
        femur_rot_back_item = R_backOpenSim @ femur_rot3_all_deform[i, :].reshape(-1, 1)
        femur_rot_back[i, :] = femur_rot_back_item.flatten()

    femurMA_rot_back = np.zeros_like(femurMA_rot3_all_deform)
    for i in range(femurMA_rot3_all_deform.shape[0]):
        femurMA_rot_back_item = R_backOpenSim @ femurMA_rot3_all_deform[i, :].reshape(-1, 1)
        femurMA_rot_back[i, :] = femurMA_rot_back_item.flatten()

    femurMarker_rot_back = np.zeros_like(femurMarker_rot3_all_deform)
    for i in range(femurMarker_rot3_all_deform.shape[0]):
        femurMarker_rot_back_item = R_backOpenSim @ femurMarker_rot3_all_deform[i, :].reshape(-1, 1)
        femurMarker_rot_back[i, :] = femurMarker_rot_back_item.flatten()

    H_ZY_back = np.dot(R_backOpenSim, H_transfer.T).T
    
    
    femur_back = np.zeros_like(femur_rot_back)
    for i in range(femur_rot_back.shape[0]):
        femur_back_item = femur_rot_back[i, :] - centroidValueLGtroch
        femur_back[i, :] = femur_back_item

    H_back = H_ZY_back + centroidValueLGtroch

    femurMA_back = np.zeros_like(femurMA_rot_back)
    for i in range(femurMA_rot_back.shape[0]):
        femurMA_back_item = femurMA_rot_back[i, :] - centroidValueLGtroch
        femurMA_back[i, :] = femurMA_back_item

    femurMarker_back = np.zeros_like(femurMarker_rot_back)
    for i in range(femurMarker_rot_back.shape[0]):
        femurMarker_back_item = femurMarker_rot_back[i, :] - centroidValueLGtroch
        femurMarker_back[i, :] = femurMarker_back_item

    femur_Rotated = femur_back
    femurMA_Rotated = femurMA_back
    femurMarker_Rotated = femurMarker_back

    fig = plt.figure(figsize=(5, 9.5), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(polys, femur_Rotated[:, 0], femur_Rotated[:, 1], femur_Rotated[:, 2],
                    edgecolor='black', linestyle=':', cmap='gray')
    ax.plot_trisurf(polys, femur_start[:, 0], femur_start[:, 1], femur_start[:, 2], edgecolor='black')
    if femurMA_Rotated.size != 0:
        ax.scatter3D(femurMA_Rotated[:, 0], femurMA_Rotated[:, 1], femurMA_Rotated[:, 2], color='red')
        ax.scatter3D(femurMuscle_start[:, 0], femurMuscle_start[:, 1], femurMuscle_start[:, 2],
                     s=20, color='blue')
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    plt.show()
    
    ##Convert to OpenSim coordinate
    femur_OpenSim = coordinatesOpenSim(femur_Rotated)
    femurMuscle_OpenSim = coordinatesOpenSim(femurMA_Rotated)
    femurMarker_OpenSim = coordinatesOpenSim(femurMarker_Rotated)
    
    
    ##change the name of the femur in the gait2392 or rajagopal model file
    
    # Try if input model == gait model (ThelenMuscles) or rajagopal (Millard - which would cause an error)
    muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle[0]  # Assuming it's a list

    if answerLeg == rightbone:
        for i in range(len(dataModel["Model"]["BodySet"]["objects"]["Body"])):
            if 'femur_r' in dataModel["Model"]["BodySet"]["objects"]["Body"][i].Attributes.name:
                if hasattr(dataModel["Model"]["BodySet"]["objects"]["Body"][i], 'attached_geometry'):
                    dataModel["Model"]["BodySet"]["objects"]["Body"][i].attached_geometry.Mesh.mesh_file = c
                else:
                    dataModel["Model"]["BodySet"]["objects"]["Body"][i].VisibleObject.GeometrySet.objects.DisplayGeometry.geometry_file.Text = c
                break
    else:
        for i in range(len(dataModel["Model"]["BodySet"]["objects"]["Body"])):
            if 'femur_l' in dataModel["Model"]["BodySet"]["objects"]["Body"][i].Attributes.name:
                if hasattr(dataModel["Model"]["BodySet"]["objects"]["Body"][i], 'attached_geometry'):
                    dataModel["Model"]["BodySet"]["objects"]["Body"][i].attached_geometry.Mesh.mesh_file = c
                else:
                    dataModel["Model"]["BodySet"]["objects"]["Body"][i].VisibleObject.GeometrySet.objects.DisplayGeometry.geometry_file.Text = c
                break
                
                
    for i in range(len(femurMuscle)):
        if len(femurPlace1[i]) == 14:
            musclenr_femur = femurNR[i]
            string_femur = femurPlace1[i]
            dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle[0].GeometryPath.PathPointSet.objects[string_femur[0:9]][0][int(string_femur[13])].location.Text = femurMuscle_OpenSim[i]
        elif len(femurPlace1[i]) == 9:
            musclenr_femur = femurNR[i]
            string_femur = femurPlace1[i]
            dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle[0].GeometryPath.PathPointSet.objects[string_femur].location.Text = femurMuscle_OpenSim[i]
        elif len(femurPlace1[i]) == 20:
            musclenr_femur = femurNR[i]
            string_femur = femurPlace1[i]
            dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle[0].GeometryPath.PathPointSet.objects[string_femur].location.Text = femurMuscle_OpenSim[i]
        else:
            musclenr_femur = femurNR[i]
            string_femur = femurPlace1[i]
            dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle[0].GeometryPath.PathPointSet.objects[string_femur[0:20]][0][int(string_femur[23])].location.Text = femurMuscle_OpenSim[i]

    for i in range(len(markerFemur_start)):
        musclenr = markerFemurNR[i]
        markerset.OpenSimDocument.MarkerSet.objects.Marker[0][musclenr].location.Text = femurMarker_OpenSim[i]


    try:
        muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle[0]
        if strcmp(answerLeg, rightbone) == 1:
            for i in range(len(dataModel["Model"]["BodySet"]["objects"]["Body"])):
                if 'femur_r' in dataModel["Model"]["BodySet"]["objects"]["Body"][i].Attributes.name:
                    if 'attached_geometry' in dataModel["Model"]["BodySet"]["objects"]["Body"][i]:
                        dataModel["Model"]["BodySet"]["objects"]["Body"][i].attached_geometry.Mesh.mesh_file = c
                    else:
                        dataModel["Model"]["BodySet"]["objects"]["Body"][i].VisibleObject.GeometrySet.objects.DisplayGeometry.geometry_file.Text = c
                    break
        else:
            for i in range(len(dataModel["Model"]["BodySet"]["objects"]["Body"])):
                if 'femur_l' in dataModel["Model"]["BodySet"]["objects"]["Body"][i].Attributes.name:
                    if 'attached_geometry' in dataModel["Model"]["BodySet"]["objects"]["Body"][i]:
                        dataModel["Model"]["BodySet"]["objects"]["Body"][i].attached_geometry.Mesh.mesh_file = c
                    else:
                        dataModel["Model"]["BodySet"]["objects"]["Body"][i].VisibleObject.GeometrySet.objects.DisplayGeometry.geometry_file.Text = c
                    break

        for i in range(len(femurMuscle)):
            if '{' in femurPlace1[i]:
                ind = femurPlace1[i].find('{')
                tmp_place = femurPlace1[i][:ind-1]
                dataModel.OpenSimDocument.Model.ForceSet.objects[femurMuscleType[i]][0].GeometryPath.PathPointSet.objects[tmp_place][0][int(femurPlace1[i][-2])].location.Text = femurMuscle_OpenSim[i]
            else:
                tmp_place = femurPlace1[i]
                dataModel.OpenSimDocument.Model.ForceSet.objects[femurMuscleType[i]][0].GeometryPath.PathPointSet.objects[tmp_place].location.Text = femurMuscle_OpenSim[i]
    except Exception as ME:
        if 'Thelen2003Muscle' in str(ME):
            # Handle the Thelen2003Muscle exception here
            pass
            
    ## change the name of the model
    type = 'deformed'
    modelNamePrint = '{}_{}'.format(modelName, type)
    dataModel.OpenSimDocument.Model.Attributes.name = 'deformed_model'  # modelNamePrint
    
    
    ## Export the whole gait2392 model file - rotated muscle attachements and correct bone rotataion names
    Model2392_rotatedfemur = struct2xml(dataModel)

    if answerLeg == rightbone:
        placeNameModel = os.path.join(direct, place, modelName + '.osim')
    else:
        placeNameModel = os.path.join(direct, place, 'FINAL_PERSONALISEDTORSIONS.osim')
        
        
    with open(placeNameModel, 'w') as f:
        f.write(Model2392_rotatedfemur)

    print('New model file has been saved')
    
    markersetup_rotatedfemur = struct2xml(markerset)

# Write the marker set as an XML file
    markersName = answerNameMarkerFemur
    markerNameOut = f"{modelName}_{markersName}"
    placeNameMarkers = f"{direct}{place}{markerNameOut}"
    with open(placeNameMarkers, 'w') as f:
        f.write(markersetup_rotatedfemur)

    print('New marker set has been saved')

    os.chdir('..')





            
    
    
    


    
    
    


 