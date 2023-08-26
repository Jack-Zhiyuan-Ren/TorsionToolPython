clear all; close all
mfile_name = mfilename('fullpath');
[pathstr,name,ext] = fileparts(mfile_name);
cd(pathstr);
try
    delete('DEFORMED_MODEL/*');
end
addpath(genpath(pwd))

%% values for Rajagopal as base model
model = 'Rajagopal/Rajagopal2015.osim';
GeometryFolder = 'Rajagopal/Geometry';
applyTibiaTorsionToJointOffset = 1;
default_Anteversion = 21; 
default_NeckShaftAngle = 121;
default_TibiaTorsion = 24;

%% Right Femur
markerset = 'MarkerSet.xml'; 

deform_bone = 'F'; 
which_leg = 'R'; 
angle_AV_right = 15.7; % right anteversion angle (in degrees) %
angle_NS_right = 126.5; % right neck-shaft angle (in degrees) %
deformed_model = ['rightNSA' num2str(angle_NS_right) '_rightAVA' num2str(angle_AV_right) ];

make_PEmodel_FemurOnly(model, deformed_model, markerset, deform_bone, which_leg, angle_AV_right - default_Anteversion, angle_NS_right - default_NeckShaftAngle, GeometryFolder);

%% left femur
model = [deformed_model '.osim']; 
markerset = [deformed_model '_' markerset]; 

deform_bone = 'F'; 
which_leg = 'L'; 
angle_AV_left = 12.5; % left anteversion angle (in degrees) %
angle_NS_left = 129.5; % left neck-shaft angle (in degrees) %
deformed_model = [ 'leftNSA' num2str(angle_NS_left) '_leftAVA' num2str(angle_AV_left)]; 
make_PEmodel_FemurOnly(model, deformed_model, markerset, deform_bone, which_leg, angle_AV_left - default_Anteversion, angle_NS_left - default_NeckShaftAngle, GeometryFolder);


%% PE Model for Femur Only
% function [ ready ] = make_PEmodel_FemurOnly( answerModel, deformed_model, answerMarkerSet, deform_bone, which_leg, angle, angle_NS, geometryFolder)
% place = [cd '\DEFORMED_MODEL\'];
% 
% % what model you want to deform
% if strcmp(which_leg, 'R') == 1 &&  strcmp(deform_bone, 'F') == 1;
%         answerModel_tmp = [ answerModel];
%         answerMarkerSet_tmp = [ answerMarkerSet];
% else
%     answerModel_tmp = [place answerModel];
%     answerMarkerSet_tmp = [place answerMarkerSet];
% end
% dataModel = xml2struct(answerModel_tmp);
% 
% if ~exist("DEFORMED_MODEL\Geometry", 'dir')
%     mkdir("DEFORMED_MODEL\Geometry")
% end
% 
% if strcmp(deform_bone, 'F') && strcmp(which_leg, 'R')
%     for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%         if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%             if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry, 'Mesh')
%                 if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                         .Mesh, 2) > 1
%                     for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                             .Mesh, 2)
%                         try
%                             vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                                 .Mesh{j}.mesh_file.Text;
%                             copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                         end
%                     end
%                 else
%                     try
%                         vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                             .Mesh.mesh_file.Text;
%                         copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                     end
%                 end
%             end
%         else
%             try
%                 if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                         .GeometrySet.objects.DisplayGeometry, 2) > 1
%                     for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                         .GeometrySet.objects.DisplayGeometry, 2)
%                         try
%                             vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                                 .GeometrySet.objects.DisplayGeometry{j}.geometry_file.Text;
%                             copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                         end
%                     end
% 
%                 else
%                     try
%                         vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                             .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
%                         copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                     end
%                 end
%             end
%         end
% 
%         if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'components')
%             if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components, 'PhysicalOffsetFrame')
%                 if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame, 'attached_geometry')
%                     if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry, 'Mesh')
%                         if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
%                                 .Mesh, 2) > 1
%                             for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
%                                     .Mesh, 2)
%                                 try
%                                     vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
%                                         .Mesh{j}.mesh_file.Text;
%                                     copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                                 end
%                             end
%                         else
%                             try
%                                 vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
%                                     .Mesh.mesh_file.Text;
%                                 copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% % what you want to name the deformed model
% answerNameModel = deformed_model;
% 
% % the marker set for this model.
% markerset = xml2struct(answerMarkerSet_tmp);
% answerLegFemur = which_leg;
% answerDegFemur = angle;
% 
% for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%     if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['femur_' lower(which_leg)])
%         if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%             femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                 .Mesh.mesh_file.Text;
%         else
%             femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                 .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
%         end
%     end
% end
% 
% if strcmp(answerLegFemur, 'R') == 1; % Rotation of the right foot
% %         FA_preAngle = 17.6;
% %         NS_preAngle = 123.3;
%         % The added anteversion angle is definded
%         angleCorrection = answerDegFemur;% - FA_preAngle;
%         FA_angle = -(angleCorrection*(pi/180));
%         NS_angle = -((angle_NS)*(pi/180));
%         % The geomerty of the right femur is imported
%         dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
%     else % Rotation of the left foot
% %         FA_preAngle = 17.6;
% %         NS_preAngle = 123.3;
%         % The added anteversion angle is definded
%         angleCorrection = answerDegFemur;% - FA_preAngle;
%         FA_angle = angleCorrection*(pi/180);
%         NS_angle = ((angle_NS)*(pi/180));
%         % The geometry of the left femur is imported
%         dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
% end
%     % the script for the rotation of the femur is called.
%     femur_ns(dataModel, markerset, answerLegFemur, 'R', FA_angle, NS_angle,...
%         answerNameModel,answerMarkerSet, dataFemur, place);
% 
% for i = 1 : numel(GeometryPathStarts)
%     section = filetext(GeometryPathStarts(i) : GeometryPathEnds(i));
%     i1 = strfind(section, '<objects');
%     i2 = strfind(section, '</objects');
%     if ~isempty(i2)
%         section = section(i1:i2);
%         nameIdx = strfind(section, 'name=');
%         absStart = 0;
%         absEnd = 0;
%         names = [];
%         for j = 1 : numel(nameIdx)
%             sep = section(nameIdx(j) + 5);
%             sep2 = strfind(section(nameIdx(j) + 6 : end), sep);
%             sep2 = sep2(1);
%             names{j} = section(nameIdx(j) + 6: nameIdx(j) + 4 + sep2);
%             tmp1 = strfind(section(nameIdx(j) - 25 : end), '<');
%             tmp1 = nameIdx(j) - 25 + tmp1(1);
%             objType = section(tmp1 : nameIdx(j)-2);
%             endIdx = strfind(section(tmp1 : end), ['</' objType]);
%             endIdx = tmp1 + endIdx;
%             section(tmp1-1 : endIdx + length(objType) +1);
%             fullText{j} = section(tmp1-1 : endIdx + length(objType) +1);
%             if j == 1
%                 absStart = tmp1-1;
%             end
%             if j == numel(nameIdx)
%                 absEnd = endIdx + length(objType) +1;
%             end
%         end
%         [~, order] = sort(names);
%         sectionNewOrder = '';
%         for j = 1 : numel(order)
%             sectionNewOrder = [sectionNewOrder fullText{order(j)}];
%         end
%         newFileText = strrep(newFileText, section(absStart : absEnd), sectionNewOrder);
%     end
% end
% 
% fid = fopen(['DEFORMED_MODEL/' answerNameModel '.osim'],'w');
% fprintf(fid, newFileText);
% fclose(fid);
% 
% end
% 
% 
% function femur_ns(dataModel, markerset, answerLeg, rightbone, FA_angle, NS_angle,answerNameModelFemur, ...
%     answerNameMarkerFemur, dataFemur, place)
% 
% %% Bone vertix
% % change the vertices into num and find the polys
% femur = str2num(dataFemur.VTKFile.PolyData.Piece.Points.DataArray.Text);
% polyText = dataFemur.VTKFile.PolyData.Piece.Polys.DataArray{1,1}.Text;
% 
% % The muscle attachments for the femur are put in one matrix.
% [femurMuscle, femurPlace1, femurNR, femurMuscleType] = femur_MA(dataModel, answerLeg, rightbone);
% %Find the markers attached to the femur in OpenSim
% [~, ~, ~, ~, markerFemur, markerFemurNR] = OpenSimMarkers(markerset, answerLeg, rightbone);
% 
% % The vertices for the bone and muscle attachements
% % are rotated to fit the coordinate system in MATLAB
% [femur_start]=coordinatesCorrection(femur);
% [femurMuscle_start] = coordinatesCorrection(femurMuscle);
% [markerFemur_start] = coordinatesCorrection(markerFemur);
% 
% wrapCnt = 0;
% wrapLocations = [];
% wrapRotations = [];
% wrapIndizes = [];
% for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%     if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1, i}.Attributes.name, ['femur_' lower(answerLeg)])
%         if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'WrapObjectSet')
%             if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet, 'objects')
%                 wrapObjectTypes = fieldnames(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects);
%                 for j = 1 : numel(wrapObjectTypes)
%                     for k = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}))
%                         wrapCnt = wrapCnt+1;
%                         wrapRotations(wrapCnt, :) = str2num(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.xyz_body_rotation.Text);
%                         wrapLocations(wrapCnt, :) = str2num(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.translation.Text);
%                         wrapIndizes(wrapCnt, :) = [i j k];
%                     end
%                 end
%             end
%         end
%         break;
%     end
% end
% 
% wrapRotations_New = wrapRotations;
% wrapLocations_New = wrapLocations;
% % wrapRotations = coordinatesCorrection(wrapRotations);
% wrapLocations = coordinatesCorrection(wrapLocations);
% 
% %The new femoral shaft axis is found as well as the inner, middle and our box to prepare for the rotation of the bone
% [innerBox, middleBox,innerBoxMA, innerBoxMarker, middleBoxMA, middleBoxMarker,femur_NewAxis, H_transfer,...
%     angleZX, angleZY, angleXY,centroidValueLGtroch, femurMA_NewAxis, femurMarker_NewAxis,...
%     Condyl_NewAxis, Shaft_proximal, Shaft_distal, CondylMA_NewAxis, ShaftMA_distal,CondylMarker_NewAxis,ShaftMarker_distal, ...
%     wrapLocations_NewAxis] = femurShaft_ns(dataModel, femur_start, answerLeg, rightbone, femurMuscle_start,markerFemur_start, wrapLocations);
% 
% %% import and convert the polies for figures
% femurSize = size(femur_start);
% polysplit = strsplit(polyText,'\n');
% poly3 = [];
% for i = 1:size(polysplit,2)
%     poly3 = [poly3; str2num(polysplit{1,i})];
% end
% polys = poly3+1;
% 
% %% Rotation; step 1
% 
% %  Here we determine where in the bone the vertix or the muslce attachments are and rotate them depending on it create the rotation matix
% Ry_NS = [cos(NS_angle) 0 sin(NS_angle); 0 1 0; -sin(NS_angle) 0 cos(NS_angle)];%neck shaft angle
% Rz_FA = [cos(FA_angle) -sin(FA_angle) 0; sin(FA_angle) cos(FA_angle) 0; 0 0 1];%femoral anteversion
% RotMatrix =  Rz_FA * Ry_NS;
% % The top of the femoral bone (head and greater trochanger) rotated around the femoral shaft axis, when the origin is in the center of the greater/lesser torchanter
% femur_rot1_all = []; innerBox_rot1 = []; polys_innerNumber = [];
% for i= 1: size(femur_NewAxis)
%     if ismember(femur_NewAxis(i,:), innerBox) == 1
%         item_innerBox =( RotMatrix * femur_NewAxis(i,:)')';
%         femur_rot1_all(i,:) = item_innerBox(:,:);
%         innerBox_rot1 = [innerBox_rot1;femur_rot1_all(i,:)];
%         polys_innerNumber = [polys_innerNumber; i];
%     else
%         femur_rot1_all(i,:) = femur_NewAxis(i,:);
%     end
% end
% %zeroMatrix will be the same size as tri3 filled with zeros
% zeroMatrix = polys * 0;
% %create a matix of ones to know when polys in shaft occur
% for ii = 1:size(polys_innerNumber,1)
%     zeroMatrix = zeroMatrix + (polys == polys_innerNumber(ii));
% end
% polys_inner = [];
% %sort out triangle that occur in the shaft
% for k = 1:size(zeroMatrix,1)
%     if sum(zeroMatrix(k,:))==3
%         polys_inner = [polys_inner; polys(k,:)]; %polys for the shaft - in the conter clockwise order
%     end
% end
% 
% femurMA_rot1_all = [];
% for i= 1: size(femurMA_NewAxis)
%     if ismember(femurMA_NewAxis(i,:), innerBoxMA) == 1
%         item_innerBoxMA =( RotMatrix * femurMA_NewAxis(i,:)')';
%         femurMA_rot1_all(i,:) = item_innerBoxMA(:,:);
%     else
%         femurMA_rot1_all(i,:) = femurMA_NewAxis(i,:);
%     end
% end
% femurMarker_rot1_all = [];
% for i= 1: size(femurMarker_NewAxis)
%     if ismember(femurMarker_NewAxis(i,:), innerBoxMarker) == 1
%         item_innerBoxMarker =( RotMatrix * femurMarker_NewAxis(i,:)')';
%         femurMarker_rot1_all(i,:) = item_innerBoxMarker(:,:);
%     else
%         femurMarker_rot1_all(i,:) = femurMarker_NewAxis(i,:);
%     end
% end
% innerBox_H_rot = (RotMatrix * H_transfer')';
% % centroidTroc_rot = (Ry_NS *Origin_transfer')';
% figure('position', [500, 50, 500, 950]); colormap([1,1,1])
% trisurf(polys, femur_rot1_all(:,1), femur_rot1_all(:,2), femur_rot1_all(:,3), 'edgecolor','black','LineStyle',':'); hold on
% trisurf(polys, femur_NewAxis(:,1), femur_NewAxis(:,2), femur_NewAxis(:,3), 'edgecolor','black');
% axis equal; set(gca,'FontSize',20); view(20,-10); xlabel('x'); ylabel('y'); zlabel('z')
% % axis equal; set(gca,'FontSize',20); view(250,40); xlabel('x'); ylabel('y'); zlabel('z')
% % scatter3(femurMarker_NewAxis(:,1),femurMarker_NewAxis(:,2), femurMarker_NewAxis(:,3),20,'b'); hold on
% % scatter3(femurMarker_rot1_all(:,1),femurMarker_rot1_all(:,2), femurMarker_rot1_all(:,3),20,'black'); hold on
% % scatter3(femurMA_NewAxis(:,1),femurMA_NewAxis(:,2), femurMA_NewAxis(:,3),20,'b'); hold on
% % scatter3(femurMA_rot1_all(:,1),femurMA_rot1_all(:,2), femurMA_rot1_all(:,3),20,'black'); hold on
% % create a new pair of axes inside current figu20
% % axes('position',[.60 .250 .35 .55])
% box on % put box around new pair of axes
% % trisurf(polys_inner,femur_NewAxis(:,1), femur_NewAxis(:,2), femur_NewAxis(:,3), 'edgecolor','black'); hold on
% % trisurf(polys_inner,femur_rot1_all(:,1), femur_rot1_all(:,2), femur_rot1_all(:,3), 'edgecolor','black','LineStyle',':'); hold on
% % scatter3(femurMarker_rot1_all(:,1),femurMarker_rot1_all(:,2), femurMarker_rot1_all(:,3),10,'r'); hold on %the femur with the new axis
% axis equal; view(30,-10); set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% % axis equal; view(250,40); set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
% 
% maxZ = max(innerBox(:, 3));
% minZ = min(innerBox(:, 3));
% 
% % Ry_FA = [cos(FA_angle) 0 sin(FA_angle);0 1 0; -sin(FA_angle) 0 cos(FA_angle)];
% % % Rx_NS = [cos(-NS_angle) -sin(-NS_angle) 0; sin(-NS_angle) cos(-NS_angle) 0; 0 0 1];%neck shaft angle
% % Rx_NS = [1 0 0; 0 cos(-NS_angle) -sin(-NS_angle); 0 sin(-NS_angle) cos(-NS_angle)];%neck shaft angle
% % 
% % RotMatrixForLocation = Rx_NS * Ry_FA;
% % 
% % % Ry_NS = [cos(NS_angle) 0 sin(NS_angle); 0 1 0; -sin(NS_angle) 0 cos(NS_angle)];%neck shaft angle
% % Rx_NS = [1 0 0; 0 cos(NS_angle) -sin(NS_angle); 0 sin(NS_angle) cos(NS_angle)];%neck shaft angle
% % % Rz_FA = [cos(FA_angle) -sin(FA_angle) 0; sin(FA_angle) cos(FA_angle) 0; 0 0 1];%femoral anteversion
% % Ry_FA = [1 0 0; 0 cos(-FA_angle) -sin(-FA_angle); 0 sin(-FA_angle) cos(-FA_angle)];%femoral anteversion
% % RotMatrixForRotation =  Ry_FA * Rx_NS;
% 
% for u = 1 : wrapCnt
%     currLoc = wrapLocations_NewAxis(u, :);
%     if currLoc(3) > minZ && currLoc(3) < maxZ % this wrapobject is in the top box --> rotate with this rotation matrix
% %         wrapLocations_New(i, :) = RotMatrixForLocation * wrapLocations_New(i, :)';
% %         wrapRotations_New(i, :) = RotMatrixForRotation * wrapRotations_New(i, :)';
%         
%         i = wrapIndizes(u, 1);
%         j = wrapIndizes(u, 2);
%         k = wrapIndizes(u, 3);
%         name = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.Attributes.name;
%         disp(['WrapObject "' name '" is in the area where the femur itself is rotated! This is not supported yet, adjust the location manually!']);
%     end
% end
% 
% %% Rotation: step 2
% % Find the highest and lowest value in this part
% high_middleBox = max(middleBox(:,3));
% low_middleBox = min(middleBox(:,3));
% low_value = zeros(1,3);
% high_value = zeros(1,3);
% for i = 1:size(middleBox,1)
%     if middleBox(i,2) <= low_middleBox
%         low_value = middleBox(i,:);
%     else
%         high_value = middleBox(i,:);
%     end
% end
% topmiddlebox = [0,0,high_middleBox];
% bottommiddlebox = [0,0,low_middleBox];
% middlevector = bottommiddlebox-topmiddlebox;
% distanceMiddleBox = norm(middlevector);
% polys_middle_Number = [];
% femur_rot2_all = []; middleBox_rot2 = []; middleBox_before_rot = [];
% for i= 1:size(femur_rot1_all)
%     if ismember(femur_rot1_all(i,:), middleBox) == 1
%         middleBox_before_rot = [middleBox_before_rot; femur_rot1_all(i,:)];
%         a_proj = [femur_rot1_all(i,1) femur_rot1_all(i,2) femur_rot1_all(i,3)] - topmiddlebox;
%         projv = dot(a_proj,middlevector)/distanceMiddleBox;
%         % The linear twist in the middle part of the bone is created
%         scalingRotVectZ = (abs(projv-distanceMiddleBox)/distanceMiddleBox)*FA_angle;
%         RotMatrix_rot2 = [cos(scalingRotVectZ) -sin(scalingRotVectZ) 0; sin(scalingRotVectZ) cos(scalingRotVectZ) 0; 0 0 1];
%         scalingRotVectY = (abs(projv-distanceMiddleBox)/distanceMiddleBox);
%         grad_rot2_middleBox  = RotMatrix_rot2*femur_rot1_all(i,:)';
%         femur_rot2_all(i,:) = grad_rot2_middleBox';
%         middleBox_rot2 = [middleBox_rot2; grad_rot2_middleBox'];
%         polys_middle_Number = [polys_middle_Number; i];
%     else
%         femur_rot2_all(i,:) = femur_rot1_all(i,:);
%     end
% end
% 
% %create a matix of ones to know when polys in shaft occur
% for ii = 1:size(polys_middle_Number,1)
%     zeroMatrix = zeroMatrix + (polys == polys_middle_Number(ii));
% end
% polys_middle = [];
% %sort out triangle that occur in the shaft
% for k = 1:size(zeroMatrix,1)
%     if sum(zeroMatrix(k,:))==3
%         polys_middle = [polys_middle; polys(k,:)]; %polys for the shaft - in the conter clockwise order
%     end
% end
% femurMA_rot2_all = []; middleBoxMA_rot2 = []; middleBoxMA_before_rot = [];
% for i= 1:size(femurMA_rot1_all)
%     if ismember(femurMA_rot1_all(i,:), middleBoxMA) == 1
%         middleBoxMA_before_rot = [middleBoxMA_before_rot; femurMA_rot1_all(i,:)];
%         a_proj = [femurMA_rot1_all(i,1) femurMA_rot1_all(i,2) femurMA_rot1_all(i,3)] - topmiddlebox;
%         projv = dot(a_proj,middlevector)/distanceMiddleBox;
%         % The linear twist in the middle part of the bone is created
%         scalingRotVectZ = (abs(projv-distanceMiddleBox)/distanceMiddleBox)*FA_angle;
%         RotMatrix_rot2 = [cos(scalingRotVectZ) -sin(scalingRotVectZ) 0; sin(scalingRotVectZ) cos(scalingRotVectZ) 0; 0 0 1];
%         scalingRotVectY = (abs(projv-distanceMiddleBox)/distanceMiddleBox);
%         grad_rot2_middleBoxMA  = RotMatrix_rot2*femurMA_rot1_all(i,:)';
%         femurMA_rot2_all(i,:) = grad_rot2_middleBoxMA';
%         middleBoxMA_rot2 = [middleBoxMA_rot2; grad_rot2_middleBoxMA'];
%     else
%         femurMA_rot2_all(i,:) = femurMA_rot1_all(i,:);
%     end
% end
% 
% femurMarker_rot2_all = []; middleBoxMarker_rot2 = []; middleBoxMarker_before_rot = [];
% for i= 1:size(femurMarker_rot1_all)
%     if ismember(femurMarker_rot1_all(i,:), middleBoxMarker) == 1
%         middleBoxMarker_before_rot = [middleBoxMarker_before_rot; femurMarker_rot1_all(i,:)];
%         a_proj = [femurMarker_rot1_all(i,1) femurMarker_rot1_all(i,2) femurMarker_rot1_all(i,3)] - topmiddlebox;
%         projv = dot(a_proj,middlevector)/distanceMiddleBox;
%         % The linear twist in the middle part of the bone is created
%         scalingRotVectZ = (abs(projv-distanceMiddleBox)/distanceMiddleBox)*FA_angle;
%         RotMatrix_rot2 = [cos(scalingRotVectZ) -sin(scalingRotVectZ) 0; sin(scalingRotVectZ) cos(scalingRotVectZ) 0; 0 0 1];
%         scalingRotVectY = (abs(projv-distanceMiddleBox)/distanceMiddleBox);
%         femurMarker_rot2_all(i,:) = grad_rot2_middleBoxMarker';
%         middleBoxMarker_rot2 = [middleBoxMarker_rot2; grad_rot2_middleBoxMarker'];
%     else
%         femurMarker_rot2_all(i,:) = femurMarker_rot1_all(i,:);
%     end
% end
% % plot the femoral bone as a scatter plot with inner and middle box rotated and twisted.
% figure('position', [1000, 50, 500, 950]); colormap([1,1,1]);
% trisurf(polys,femur_rot2_all(:,1),femur_rot2_all(:,2),femur_rot2_all(:,3), 'edgecolor','black','LineStyle',':'); hold on
% trisurf(polys,femur_NewAxis(:,1),femur_NewAxis(:,2),femur_NewAxis(:,3), 'edgecolor','black');
% % scatter3(femurMA_NewAxis(:,1),femurMA_NewAxis(:,2), femurMA_NewAxis(:,3),20,'b'); hold on
% % scatter3(femurMA_rot2_all(:,1),femurMA_rot2_all(:,2), femurMA_rot2_all(:,3),20,'black'); hold on
% % scatter3(femurMarker_NewAxis(:,1),femurMarker_NewAxis(:,2), femurMarker_NewAxis(:,3),20,'b'); hold on
% % scatter3(femurMarker_rot2_all(:,1),femurMarker_rot2_all(:,2), femurMarker_rot2_all(:,3),20,'black'); hold on
% axis equal; set(gca,'FontSize',20); view(30,-10); grid on; xlabel('x'); ylabel('y'); zlabel('z');
% % axis equal; set(gca,'FontSize',20); view(250,40); grid on; xlabel('x'); ylabel('y'); zlabel('z');
% set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
% set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', [])
% % create a new pair of axes inside current figure
% axes('position',[.65 .250 .35 .55])
% box on % put box around new pair of axes
% trisurf(polys_middle,femur_rot2_all(:,1), femur_rot2_all(:,2), femur_rot2_all(:,3), 'edgecolor','black','LineStyle',':'); hold on
% trisurf(polys_middle,femur_NewAxis(:,1), femur_NewAxis(:,2), femur_NewAxis(:,3), 'edgecolor','black'); hold on
% trisurf(polys_inner,femur_NewAxis(:,1), femur_NewAxis(:,2), femur_NewAxis(:,3), 'edgecolor','black'); hold on
% trisurf(polys_inner,femur_rot2_all(:,1), femur_rot2_all(:,2), femur_rot2_all(:,3), 'edgecolor','black','LineStyle',':'); hold on
% axis equal; view(30,-10); set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% % axis equal; view(250,40); set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
% 
% 
% maxZ = max(middleBox(:, 3));
% minZ = min(middleBox(:, 3));
% for u = 1 : wrapCnt
%     currLoc = wrapLocations_NewAxis(u, :);
%     if currLoc(3) > minZ && currLoc(3) < maxZ % this wrapobject is in the middle box --> rotate with this rotation matrix 2
%         %         wrapRotations_New(i, :) = RotMatrix_rot2 * wrapRotations_New(i, :)';
%         %         wrapLocations_New(i, :) = RotMatrix_rot2 * wrapLocations_New(i, :)';
%         i = wrapIndizes(u, 1);
%         j = wrapIndizes(u, 2);
%         k = wrapIndizes(u, 3);
%         name = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.Attributes.name;
%         disp(['WrapObject "' name '" is in the area where the femur itself is rotated! This is not supported yet, adjust the location manually!']);
%     end
% end
% 
% 
% %% Rotation; Step 3
% %MOVE THE FEMORAL HEAD BACK TO FIT CONDYLAR WITHOUT MOVING THE CONDYLAR
% %start by transfering the bone back to the H_ZX
% transfer_step3 = H_transfer-innerBox_H_rot;
% 
% %find the distance between the highest value in the distal shaft to the highese value in the condylar
% maxCondyl = max(Condyl_NewAxis(:,3));
% condyl_top = zeros(1,3);
% for i = 1:size(Condyl_NewAxis,1)
%     if Condyl_NewAxis(i,3) >= maxCondyl
%         condyl_top = Condyl_NewAxis(i,:);
%     end
% end
% 
% minShaft_prox = min(Shaft_proximal(:,3));
% shaft_prox_min = zeros(1,3);
% for i = 1:size(Shaft_proximal,1)
%     if Shaft_proximal(i,3) <= minShaft_prox
%         shaft_prox_min = Shaft_proximal(i,:);
%     end
% end
% distance_shaft_distal = norm( condyl_top- shaft_prox_min);
% 
% %
% femur_rot3_all = []; % The outer box is moved back to restore the position of the femoral head to the acetabulum
% femur_rot3_all_deform = []; % the oter box has been restored and the distal part of the femoral shaft is gradually deformed to fit the condylar (which does not move)
% for i = 1:size(femur_rot2_all,1)
%     if ismember(femur_rot2_all(i,:), Condyl_NewAxis) == 1 % the condylar do not move but the rest is translated restore the femoral head
%         femur_rot3_all(i,:) = femur_rot2_all(i,:);
%     else
%         item_rot3 =femur_rot2_all(i,:) + transfer_step3; %The outer box are transfered back to the postion of the femoral head
%         femur_rot3_all(i,:) = item_rot3(:,:);
%     end
%     if ismember(femur_rot2_all(i,:), Shaft_distal) == 1
%         scaler = (abs(norm(femur_rot3_all(i,:)-condyl_top))-distance_shaft_distal)/(distance_shaft_distal) * transfer_step3;
%         item_rot3_distal = femur_rot3_all(i,:)+ scaler;
%         femur_rot3_all_deform(i,:) = item_rot3_distal(:,:);
%     else
%         femur_rot3_all_deform(i,:) = femur_rot3_all(i,:);
%     end
% end
% % Treat the muscles attachements in the same way as the bone verticies.
% % centroidTroc_rot3 = centroidTroc_rot + transfer_step3;
% femurMA_rot3_all = []; femurMA_rot3_all_deform = []; test = [];
% for i = 1:size(femurMA_rot2_all,1)
%     if ismember(femurMA_rot2_all(i,:), CondylMA_NewAxis) == 1 % the condylar do not move
%         femurMA_rot3_all(i,:) = femurMA_rot2_all(i,:);
%         test = [test; femurMA_rot2_all(i,:)];
%     else
%         item_rot3 =femurMA_rot2_all(i,:) + transfer_step3; %The inner and middel box are transfered back to the postion of the femoral head
%         femurMA_rot3_all(i,:) = item_rot3(:,:);
%     end
%     if ismember(femurMA_rot2_all(i,:), ShaftMA_distal) == 1
%         scaler = (abs(norm(femurMA_rot3_all(i,:)-condyl_top))-distance_shaft_distal)/(distance_shaft_distal) * transfer_step3;
%         item_rot3_distal = femurMA_rot3_all(i,:)+ scaler;
%         femurMA_rot3_all_deform(i,:) = item_rot3_distal(:,:);
%     else
%         femurMA_rot3_all_deform(i,:) = femurMA_rot3_all(i,:);
%     end
% end
% % treat the markers in the same way as the bone verticies
% femurMarker_rot3_all = []; femurMarker_rot3_all_deform = [];
% for i = 1:size(femurMarker_rot2_all,1)
%     if ismember(femurMarker_rot2_all(i,:), CondylMarker_NewAxis) == 1 % the condylar do not move
%         femurMarker_rot3_all(i,:) = femurMarker_rot2_all(i,:);
%     else
%         item_rot3 =femurMarker_rot2_all(i,:) + transfer_step3; %The inner and middel box are transfered back to the postion of the femoral head
%         femurMarker_rot3_all(i,:) = item_rot3(:,:);
%     end
%     if ismember(femurMarker_rot2_all(i,:), ShaftMarker_distal) == 1
%         scaler = (abs(norm(femurMarker_rot3_all(i,:)-condyl_top))-distance_shaft_distal)/(distance_shaft_distal) * transfer_step3;
%         item_rot3_distal = femurMarker_rot3_all(i,:)+ scaler;
%         femurMarker_rot3_all_deform(i,:) = item_rot3_distal(:,:);
%     else
%         femurMarker_rot3_all_deform(i,:) = femurMarker_rot3_all(i,:);
%     end
% end
% 
% maxShaft_prox = max(Shaft_proximal(:,3));
% % treat the wrapping objects in the same way
% for i = 1 : wrapCnt
%     currLoc = wrapLocations_NewAxis(i, :);
%     wrapRotations_New(i, :) = wrapRotations_New(i, :);
%     if currLoc(3) < maxCondyl % it is at the condyles, don't do any rotation
% %         wrapRotations_New(i, :) = wrapRotations_New(i, :);
%         wrapLocations_New(i, :) = wrapLocations_New(i, :);
%     else
% %         wrapRotations_New(i, :) = wrapRotations_New(i, :) + transfer_step3;
% %         wrapLocations_New(i, :) = wrapLocations_New(i, :) + transfer_step3;
% %         wrapRotations_New(i, :) = wrapRotations_New(i, :) + transfer_step3;
%         wrapLocations_New(i, :) = wrapLocations_New(i, :) + coordinatesOpenSim(transfer_step3);
%     end
%     if currLoc(3) < maxShaft_prox && currLoc(3) > minShaft_prox
%         scaler = (abs(norm(wrapLocations_New(i,:)-condyl_top))-distance_shaft_distal)/(distance_shaft_distal) * transfer_step3;
% %         wrapRotations_New(i, :) = wrapRotations_New(i, :) + scaler;
%         wrapLocations_New(i, :) = wrapLocations_New(i, :) + coordinatesOpenSim(scaler);
%     end
% end
% 
% % wrapRotations_New = coordinatesOpenSim(wrapRotations_New);
% % wrapLocations_New = coordinatesOpenSim(wrapLocations_New);
% for u = 1 : wrapCnt
%     i = wrapIndizes(u, 1);
%     j = wrapIndizes(u, 2);
%     k = wrapIndizes(u, 3);
%     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.xyz_body_rotation.Text = num2str(wrapRotations_New(u, :));
%     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.WrapObjectSet.objects.(wrapObjectTypes{j}){1, k}.translation.Text = num2str(wrapLocations_New(u, :));
% end
% 
% % plot the femoral bone as a scatter plot with inner and middle box rotated and twisted.
% % figure('position', [1000, 50, 500, 950]); colormap([1,1,1]);
% % trisurf(polys,femur_rot3_all_deform(:,1),femur_rot3_all_deform(:,2),femur_rot3_all_deform(:,3), 'edgecolor','black','LineStyle',':'); hold on
% % trisurf(polys,femur_NewAxis(:,1),femur_NewAxis(:,2),femur_NewAxis(:,3), 'edgecolor','black');
% % % scatter3(Shaft_distal(:,1),Shaft_distal(:,2), Shaft_distal(:,3),10,'r'); hold on %the femur with the new axis
% if ~isempty(femurMA_NewAxis)
%     scatter3(femurMA_NewAxis(:,1),femurMA_NewAxis(:,2), femurMA_NewAxis(:,3),20,'b'); hold on
%     scatter3(femurMA_rot3_all_deform(:,1),femurMA_rot3_all_deform(:,2), femurMA_rot3_all_deform(:,3),20,'black'); hold on
% end
% % scatter3(femurMarker_NewAxis(:,1),femurMarker_NewAxis(:,2), femurMarker_NewAxis(:,3),20,'b'); hold on
% % scatter3(femurMarker_rot3_all_deform(:,1),femurMarker_rot3_all_deform(:,2), femurMarker_rot3_all_deform(:,3),20,'black'); hold on
% % axis equal; set(gca,'FontSize',16); view(30,-10); grid on; xlabel('x'); ylabel('y'); zlabel('z')
% 
% angleZX_back = - angleZX;
% angleZY_back = -angleZY;
% 
% Rx_FA = [ 1 0 0;0 cos(angleZY_back) -sin(angleZY_back); 0 sin(angleZY_back) cos(angleZY_back)];
% if strcmp(answerLeg, rightbone) == 1;
%     Ry_FA = [cos(angleZX_back)  0 sin(angleZX_back); 0 1 0; -sin(angleZX_back) 0 cos(angleZX_back)];
% else
%     Ry_FA = [cos(-angleZX_back)  0 sin(-angleZX_back); 0 1 0; -sin(-angleZX_back) 0 cos(-angleZX_back)];
% end
% 
% angleXY_back = -angleXY;
% if strcmp(answerLeg, rightbone) == 1;
%     Rz_FA=[cos(angleXY_back) -sin(angleXY_back) 0; sin(angleXY_back) cos(angleXY_back) 0; 0 0 1];
% else
%     Rz_FA=[cos(-angleXY_back) -sin(-angleXY_back) 0; sin(-angleXY_back) cos(-angleXY_back) 0; 0 0 1];
% end
% R_backOpenSim = Rx_FA*Ry_FA*Rz_FA;
% 
% femur_rot_back = zeros(size(femur_rot3_all_deform,1),3);
% for i = 1:size(femur_rot3_all_deform,1)
%     femur_rot_back_item = R_backOpenSim * femur_rot3_all_deform(i,:)';
%     femur_rot_back(i,:) = femur_rot_back_item';
% end
% femurMA_rot_back = zeros(size(femurMA_rot3_all_deform,1),3);
% for i = 1:size(femurMA_rot3_all_deform,1)
%     femurMA_rot_back_item = R_backOpenSim * femurMA_rot3_all_deform(i,:)';
%     femurMA_rot_back(i,:) = femurMA_rot_back_item';
% end
% femurMarker_rot_back = zeros(size(femurMarker_rot3_all_deform,1),3);
% for i = 1:size(femurMarker_rot3_all_deform,1)
%     femurMarker_rot_back_item = R_backOpenSim * femurMarker_rot3_all_deform(i,:)';
%     femurMarker_rot_back(i,:) = femurMarker_rot_back_item';
% end
% H_ZY_back = (R_backOpenSim * H_transfer')';
% % centroidTroc_ZY_back = (R_backOpenSim *centroidTroc_rot3')';
% 
% %%
% femur_back = [];
% for i = 1:size(femur_rot_back,1)
%     femur_back_item = femur_rot_back(i,:)- centroidValueLGtroch;
%     femur_back = [femur_back; femur_back_item];
% end
% H_back = H_ZY_back + centroidValueLGtroch;
% % centroidTroc_back = centroidTroc_ZY_back + centroidValueLGtroch;
% femurMA_back = [];
% for i = 1:size(femurMA_rot_back,1)
%     femurMA_back_item = femurMA_rot_back(i,:)- centroidValueLGtroch;
%     femurMA_back = [femurMA_back; femurMA_back_item];
% end
% femurMarker_back = [];
% for i = 1:size(femurMarker_rot_back,1)
%     femurMarker_back_item = femurMarker_rot_back(i,:)- centroidValueLGtroch;
%     femurMarker_back = [femurMarker_back; femurMarker_back_item];
% end
% femur_Rotated = femur_back;
% femurMA_Rotated = femurMA_back;
% femurMarker_Rotated = femurMarker_back;
% 
% figure('position', [1400, 50, 500, 950])
% colormap([1,1,1]);
% trisurf(polys,femur_Rotated(:,1),femur_Rotated(:,2), femur_Rotated(:,3), 'edgecolor','black','LineStyle',':'); hold on
% trisurf(polys,femur_start(:,1),femur_start(:,2), femur_start(:,3), 'edgecolor','black');
% if ~isempty(femurMA_Rotated)
%     scatter3(femurMA_Rotated(:,1),femurMA_Rotated(:,2), femurMA_Rotated(:,3), 'red');
%     scatter3(femurMuscle_start(:,1),femurMuscle_start(:,2), femurMuscle_start(:,3),20,'b'); hold on
% end
% 
% % scatter3(femurMarker_back(:,1),femurMarker_back(:,2), femurMarker_back(:,3), 'red')
% % grid on; axis equal; set(gca,'FontSize',20); view(250,40); xlabel('x'); ylabel('y'); zlabel('z');
% grid on; axis equal; set(gca,'FontSize',20); view(30,-10); xlabel('x'); ylabel('y'); zlabel('z');
% set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
% set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
% 
% %% opensim coordinates
% femur_OpenSim = coordinatesOpenSim(femur_Rotated);
% femurMuscle_OpenSim = coordinatesOpenSim(femurMA_Rotated);
% femurMarker_OpenSim = coordinatesOpenSim(femurMarker_Rotated);
% 
% %% export the files again as xml files
% % convert the femur data back to string
% femur_rotated =sprintf('\t\t\t%+8.6f %+8.6f %+8.6f\n',femur_OpenSim');
% femur_Repared = strrep(femur_rotated,'+',' ');
% % replace the generic data with the rotated bone
% dataFemur.VTKFile.PolyData.Piece.Points.DataArray.Text = femur_Repared;
% 
% % convert the struct back to xml file
% Femur_rotated = struct2xml(dataFemur);
% %name and placement of the femoral bone file
% direct = [];
% % export - write the model as an xml  - remember to save as a vtp file
% if strcmp(answerLeg, rightbone) == 1;
%     modelName = answerNameModelFemur;
%     boneName = 'femurR_rotated.vtp';
%     c = sprintf('%s_%s' ,modelName,boneName);
%     placeNameFemur = sprintf('%s',direct, place, c);
%     %write the model as an xml file
%     FID_femurR = fopen(placeNameFemur,'w');
%     fprintf(FID_femurR,Femur_rotated);
%     fclose(FID_femurR);
% else
%     modelName = answerNameModelFemur;
%     boneName = 'femurL_rotated.vtp';
%     c = sprintf('%s_%s' ,modelName,boneName);
%     placeNameFemur = sprintf('%s',direct, place, c);
%     FID_femurL = fopen(placeNameFemur,'w');
%     fprintf(FID_femurL,Femur_rotated);
%     fclose(FID_femurL);
% end
% 
% %% change the name of the femur in the gait2392 or rajagopal model file
% try
%     %try if input model == gait model (ThelenMuscles) or
%     %rajagopal (Millard - which would cause an error)
%     muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,1};
%     if strcmp(answerLeg, rightbone) == 1;
%         for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%             if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1, i}.Attributes.name, 'femur_r')
%                 if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%                     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                         .Mesh.mesh_file = c;
%                 else
%                     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                         .GeometrySet.objects.DisplayGeometry.geometry_file.Text = c;
%                 end
%                 break;
%             end
%         end
%     else
%         for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%             if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1, i}.Attributes.name, 'femur_l')
%                 if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%                     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                         .Mesh.mesh_file = c;
%                 else
%                     dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                         .GeometrySet.objects.DisplayGeometry.geometry_file.Text = c;
%                 end
%                 break;
%             end
%         end
%     end
%     for i = 1:size(femurMuscle,1)
%         if size(femurPlace1{i,1},2) == 14 ;
%             musclenr_femur = femurNR(i,:);
%             string_femur = femurPlace1{i,:};
%             dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,musclenr_femur}...
%                 .GeometryPath.PathPointSet.objects.(string_femur(1:9)){1,str2num(string_femur(13))}.location.Text = femurMuscle_OpenSim(i,:);
%         elseif size(femurPlace1{i,1},2) == 9;
%             musclenr_femur = femurNR(i,:);
%             string_femur = femurPlace1{i,:};
%             dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,musclenr_femur}...
%                 .GeometryPath.PathPointSet.objects.(string_femur).location.Text = femurMuscle_OpenSim(i,:);
%         elseif size(femurPlace1{i,1},2) == 20;
%             musclenr_femur = femurNR(i,:);
%             string_femur = femurPlace1{i,:};
%             dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,musclenr_femur}...
%                 .GeometryPath.PathPointSet.objects.(string_femur).location.Text = femurMuscle_OpenSim(i,:);
%         else
%             musclenr_femur = femurNR(i,:);
%             string_femur = femurPlace1{i,:};
%             dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,musclenr_femur}...
%                 .GeometryPath.PathPointSet.objects.(string_femur(1:20)){1,str2num(string_femur(24))}.location.Text = femurMuscle_OpenSim(i,:);
%         end
%     end
%     for i = 1:size(markerFemur_start,1)
%         musclenr = markerFemurNR(i,:);
%         markerset.OpenSimDocument.MarkerSet.objects.Marker{1,musclenr}.location.Text = femurMarker_OpenSim(i,:);
%     end
% catch ME
%     if contains(ME.message,'Thelen2003Muscle')
%         muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,1};
%         if strcmp(answerLeg, rightbone) == 1;
%             for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%                 if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1, i}.Attributes.name, 'femur_r')
%                     if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%                         dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                             .Mesh.mesh_file = c;
%                     else
%                         dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                             .GeometrySet.objects.DisplayGeometry.geometry_file.Text = c;
%                     end
%                     break;
%                 end
%             end
%         else
%             for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
%                 if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1, i}.Attributes.name, 'femur_l')
%                     if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
%                         dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
%                             .Mesh.mesh_file = c;
%                     else
%                         dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
%                             .GeometrySet.objects.DisplayGeometry.geometry_file.Text = c;
%                     end
%                     break;
%                 end
%             end
%         end
% 
% 
%         for i = 1:size(femurMuscle,1)
%             if contains(femurPlace1{i}, '{')
%                 ind = strfind(femurPlace1{i}, '{');
%                 tmp_place = femurPlace1{i}(1 : ind-1);
%                 dataModel.OpenSimDocument.Model.ForceSet.objects.(femurMuscleType{i}){1, femurNR(i)}...
%                     .GeometryPath.PathPointSet.objects.(tmp_place){1,str2double(femurPlace1{i}(end-1))}.location.Text = femurMuscle_OpenSim(i,:);
%             else
%                 tmp_place = femurPlace1{i};
%                 dataModel.OpenSimDocument.Model.ForceSet.objects.(femurMuscleType{i}){1, femurNR(i)}...
%                     .GeometryPath.PathPointSet.objects.(tmp_place).location.Text = femurMuscle_OpenSim(i,:);
%             end
% 
% %             dataModel.OpenSimDocument.Model.ForceSet.objects.(femurMuscleType{i}){1, femurNR(i)}...
% %                 .GeometryPath.PathPointSet.objects.(string_femur(1:9)){1,str2double(femurPlace1{i}(end-1))}.location.Text = femurMuscle_OpenSim(i,:);
% % 
% %             if size(femurPlace1{i,1},2) == 14 ; % string equals PathPoint{x,x}
% %                 musclenr_femur = femurNR(i,:);
% %                 string_femur = femurPlace1{i,:};
% %                 dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,musclenr_femur}...
% %                     .GeometryPath.PathPointSet.objects.(string_femur(1:9)){1,str2num(string_femur(13))}.location.Text = femurMuscle_OpenSim(i,:);
% %             elseif size(femurPlace1{i,1},2) == 9;
% %                 musclenr_femur = femurNR(i,:);
% %                 string_femur = femurPlace1{i,:};
% % 
% %                 dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,musclenr_femur}...
% %                     .GeometryPath.PathPointSet.objects.(string_femur).location.Text = femurMuscle_OpenSim(i,:);
% % 
% % 
% %             elseif size(femurPlace1{i,1},2) == 20;
% %                 musclenr_femur = femurNR(i,:);
% %                 string_femur = femurPlace1{i,:};
% %                 dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,musclenr_femur}...
% %                     .GeometryPath.PathPointSet.objects.(string_femur).location.Text = femurMuscle_OpenSim(i,:);
% %             else
% %                 musclenr_femur = femurNR(i,:);
% %                 string_femur = femurPlace1{i,:};
% %                 dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,musclenr_femur}...
% %                     .GeometryPath.PathPointSet.objects.(string_femur(1:20)){1,str2num(string_femur(24))}.location.Text = femurMuscle_OpenSim(i,:);
% %             end
%         end
%         for i = 1:size(markerFemur_start,1)
%             musclenr = markerFemurNR(i,:);
%             markerset.OpenSimDocument.MarkerSet.objects.Marker{1,musclenr}.location.Text = femurMarker_OpenSim(i,:);
%         end
%     else
%         error('Unknown Error in femur_ns regarding ForceSet-Muscles to identify model basis (gait or rajagopal?) check try-catch blocks')
%     end
% end
% 
% %% change the name of the model
% type= 'deformed';
% modelNamePrint = sprintf('%s_%s' ,modelName,type);
% dataModel.OpenSimDocument.Model.Attributes.name = 'deformed_model';%modelNamePrint;
% 
% %% Export the whole gait2392 model file - rotated muscle attachements and correct bone rotataion names
% % export the gait2392
% cd functions
% Model2392_rotatedfemur = struct2xml(dataModel);
% %name and placement of the femoral bone file
% placeNameModel = sprintf('%s',direct, place, modelName, '.osim');
% %write the model as an xml file
% FID_model = fopen(placeNameModel,'w');
% fprintf(FID_model,Model2392_rotatedfemur);
% fclose(FID_model);
% 
% disp('New model file has been saved')
% 
% % export the the marker setup for the scaling tool in opensim
% markersetup_rotatedfemur = struct2xml(markerset);
% %write the model as an xml file
% markersName= answerNameMarkerFemur;
% markerNameOut = sprintf('%s_%s' ,modelName,markersName);
% %name and placement of the femoral bone file
% placeNameMarkers = sprintf('%s', direct, place, markerNameOut);
% FID_markers = fopen(placeNameMarkers,'w');
% fprintf(FID_markers,markersetup_rotatedfemur);
% fclose(FID_markers);
% disp('New marker set has been saved')
% 
% cd ..
% 
% end
% 
% function [femurMuscle, femurPlace1, femurNR, femurMuscleType] = femur_MA(dataModel, answerLeg, rightbone)
% 
% %% Find the muscles in the model
% muscleType = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);
% muscleType = muscleType{1};
% muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);
% 
% %Identify the left and right leg
% if strcmp(answerLeg, rightbone) == 1;
%     femurMA = 'femur_r';
% else
%     femurMA = 'femur_l';
% end
% 
% femurMuscle =[];
% femurPlace1 = {};
% femurNR = [];
% femurMuscleType = {};
% 
% muscleTypes = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);
% 
% for k = 1 : numel(muscleTypes)
%     muscleType = muscleTypes{k};
%     muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);
% 
%     for i = 1:size(muscles,2)
%         if isfield(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}, 'GeometryPath')
%             AttachmentSize = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}.GeometryPath.PathPointSet.objects);
%         else
%             AttachmentSize = [];
%         end
%         
%         for j = 1 : numel(AttachmentSize)
%             MuscleAttachments = muscles{1,i}.GeometryPath.PathPointSet.objects.(AttachmentSize{j});
%     
%             if size(MuscleAttachments, 2) == 1 % only one PathPoint of this kind in this muscle (e.g. rect_fem PathPoint in gait2392)
%                 try
%                     CompareStrings1_femur = strcmp(femurMA, MuscleAttachments.body.Text);
%                 catch ME
%                     if contains(ME.message, 'body')
%                         CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments.socket_parent_frame.Text,'/bodyset/',''));
%                     else
%                         error('Unknown Error in femur_MA check try-catch blocks')
%                     end
%                 end
%                 if CompareStrings1_femur == 1;
%                     try
%                         femurMuscle = [femurMuscle; str2num(MuscleAttachments.location.Text)];
%                     catch
%                         femurMuscle = [femurMuscle; str2num(MuscleAttachments.socket_parent_frame.Text)];
%                     end
%                     femurMuscleType = [femurMuscleType; muscleType];
%                     femurNR = [femurNR; i];
%                     place1 = AttachmentSize{j};
%                     femurPlace1 = [femurPlace1;place1];
%                 end
%             else    % more of a kind - we need to add {1,%d} after the type
%                 for ii = 1:size(MuscleAttachments,2)
%                     try
%                         CompareStrings1_femur = strcmp(femurMA, MuscleAttachments{1,ii}.body.Text);
%                     catch ME
%                         if contains(ME.message, 'body')
%                             CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments{1,ii}.socket_parent_frame.Text,'/bodyset/',''));
%                         else
%                             error('Unknown Error in femur_MA check try-catch blocks')
%                         end
%                     end
%                     if CompareStrings1_femur == 1;
%                         try
%                             femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.location.Text)];
%                         catch
%                             femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.socket_parent_frame.Text)];
%                         end
%                         femurMuscleType = [femurMuscleType; muscleType];
%                         femurNR = [femurNR; i];
%                         place1 = [AttachmentSize{j} '{1,' num2str(ii) '}'];
%                         femurPlace1 = [femurPlace1;place1];
%                     end
%                 end
%             end
%         end
% 
% 
% 
%     end
% end
% 
% end
% 
% function [innerBox, middleBox, innerBoxMA, innerBoxMarker, middleBoxMA, middleBoxMarker,femurShaftLocRot, headShaftRot,...
%     angleZX, angleZY,angleXY,translationDis, femurShaftLocRotMA, femurShaftLocRotMarkers,...
%     Condylar,ShaftProx,ShaftDist, CondylarMA,ShaftDistMA,CondylarMarkers, ShaftMarkers, wrapLocationsRot]...
%     = femurShaft_ns(dataModel, femur_start, answerLeg, rightbone, femurMuscle_start,markerFemur_start, wrapLocations)
% 
% % points defining the axes, as determined by Lorenzo's code using the undeformed generic bones
% if strcmp(answerLeg, rightbone) == 1; % right leg
%     SEL = [-28.5309 5.3055 -3.3018]./1000; % saddle point femoral neck
%     SEL_epi = [-1.799 -19.7590 -418.0754]./1000; % saddle point between the two epicondyles
%     HC = [-0.1583 -0.2439 0.0038]./1000; % centre of femoral head
%     ISTHMUS = [-17.2534 4.2462 -14.4892]./1000; % centre of the femoral neck
% else % left leg
%     % at least make it symmetrical... 
%     SEL = [28.5309 5.3055 -3.3018]./1000; % saddle point femoral neck
%     SEL_epi = [1.799 -19.7590 -418.0754]./1000; % saddle point between the two epicondyles
%     HC = [0.1583 -0.2439 0.0038]./1000; % centre of femoral head
%     ISTHMUS = [17.2534 4.2462 -14.4892]./1000; % centre of the femoral neck
% 
% %     SEL = [28.5518 5.2971 -3.2637]./1000;
% %     SEL_epi = [1.8459 -19.72 -418.0768]./1000;
% %     HC = [0.1536 -0.2856 0.0421]./1000;
% %     ISTHMUS = [17.2422 2.2442 -14.4670]./1000;
% end
% 
% point1=SEL;
% point2=SEL_epi;
% t=0:.001:1;
% C=repmat(point1,length(t),1)'+(point2-point1)'*t;
% SEL_point =C(:,111)'; % select rotation point on shaft axis, at height of smaller trochanter
% % figure
% % scatter3(femur_start(:,1),femur_start(:,2),femur_start(:,3),'black'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % scatter3(0,0,0,'red');
% % scatter3(HC(1),HC(2),HC(3),'blue')
% % scatter3(ISTHMUS(1),ISTHMUS(2),ISTHMUS(3),'blue')
% % scatter3(SEL_point(1),SEL_point(2),SEL_point(3),'green')
% % scatter3(SEL(1),SEL(2),SEL(3),'blue')
% % scatter3(SEL_epi(1),SEL_epi(2),SEL_epi(3),'blue')
% % plot3([SEL_epi(1),SEL(1)],[SEL_epi(2),SEL(2)],[SEL_epi(3),SEL(3)], 'blue', 'Linewidth',3) % shaft axis
% % plot3([ISTHMUS(1),HC(1)],[ISTHMUS(2),HC(2)],[ISTHMUS(3),HC(3)], 'blue', 'Linewidth',3) % neck axis
% 
% 
% %% Transform the new femoral shaft axis to the rotation point (SEL_point)
% % move the zero point
% translationDis = [0 0 0] - SEL_point;
% %the bone
% femurShaftLoc = [];
% for i = 1:size(femur_start,1)
%     item  = femur_start(i,:)+translationDis;
%     femurShaftLoc = [femurShaftLoc; item];
% end
% 
% SEL_epiShaft = SEL_epi + translationDis;
% SEL_pointShaft = SEL_point + translationDis;
% headShaft = HC +translationDis;
% isthmusShaft = ISTHMUS + translationDis;
% 
% % the muscle attachments
% femurShaftLocMA = [];
% for i = 1:size(femurMuscle_start,1)
%     item  = femurMuscle_start(i,:)+translationDis;
%     femurShaftLocMA = [femurShaftLocMA; item];
% end
% % the markers
% femurShaftLocMarkers = [];
% for i = 1:size(markerFemur_start,1)
%     item  = markerFemur_start(i,:)+translationDis;
%     femurShaftLocMarkers = [femurShaftLocMarkers; item];
% end
% 
% % figure
% % % scatter3(femur_start(:,1),femur_start(:,2),femur_start(:,3),'black'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % scatter3(femurShaftLoc(:,1),femurShaftLoc(:,2),femurShaftLoc(:,3),'red'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % scatter3(SEL_epiShaft(:,1),SEL_epiShaft(:,2),SEL_epiShaft(:,3),'filled');
% % scatter3(SEL_pointShaft(:,1),SEL_pointShaft(:,2),SEL_pointShaft(:,3),'filled');
% % scatter3(headShaft(:,1),headShaft(:,2),headShaft(:,3),'filled');
% % scatter3(isthmusShaft(:,1),isthmusShaft(:,2),isthmusShaft(:,3),'filled');
% % plot3([SEL_epiShaft(:,1),SEL_pointShaft(:,1)],[SEL_epiShaft(:,2),SEL_pointShaft(:,2)],[SEL_epiShaft(:,3),SEL_pointShaft(:,3)], 'black', 'Linewidth',1.5) % shaft axis
% % plot3([headShaft(:,1),isthmusShaft(:,1)],[headShaft(:,2),isthmusShaft(:,2)],[headShaft(:,3),isthmusShaft(:,3)], 'black', 'Linewidth',1.5) % neck axis
% % plot3([0,0],[0,0],[-0.4,0], 'black', 'Linewidth',1.5) % z-axis knee joint
% 
% 
% aZY = [SEL_epiShaft(1,2), SEL_epiShaft(1,3)] -[0, 0];
% bZY = [0, -0.4]-[0, 0];
% angleZY = (acos(dot(aZY, bZY)/(norm(aZY)*norm(bZY))));
% Rx = [1 0 0; 0 cos(angleZY) -sin(angleZY); 0 sin(angleZY) cos(angleZY)];
% 
% aZX = [SEL_epiShaft(1,1), SEL_epiShaft(1,3)]-[0, 0];
% bZX = [0, 0.4]-[0,0];
% angleZX = pi - (acos(dot(aZX, bZX)/(norm(aZX)*norm(bZX))));
% 
% % the rotation matrix around the y -axis
% if strcmp(answerLeg, rightbone) == 1;
%     Ry = [cos(angleZX)  0 sin(angleZX); 0 1 0; -sin(angleZX) 0 cos(angleZX)];
% else
%     %This is for the left femur
%     Ry= [cos(-angleZX)  0 sin(-angleZX); 0 1 0; -sin(-angleZX) 0 cos(-angleZX)];
% end
% 
% R_transfer = Ry*Rx;
% 
% % align neck axis with x-axis
% tmp=(R_transfer*(headShaft-isthmusShaft)')';
% aXY = [tmp(1) tmp(2)];
% 
% if strcmp(answerLeg, rightbone) == 1; % right leg
%     bXY = [0.4 0];
%     angleXY=(acos(dot(aXY, bXY)/(norm(aXY)*norm(bXY))));
%     Rz=[cos(angleXY) -sin(angleXY) 0; sin(angleXY) cos(angleXY) 0; 0 0 1];
%     
% else
%     bXY = [-0.4 0];
%     angleXY=(acos(dot(aXY, bXY)/(norm(aXY)*norm(bXY))));
%     Rz=[cos(-angleXY) -sin(-angleXY) 0; sin(-angleXY) cos(-angleXY) 0; 0 0 1];
% end
% 
% R_transfer = Rz*Ry*Rx;
% 
% % the bone
% femurShaftLocRot = [];
% for i = 1:size(femur_start,1)
%     item  = (R_transfer* femurShaftLoc(i,:)')';
%     femurShaftLocRot = [femurShaftLocRot; item];
% end
% 
% wrapLocationsRot = [];
% for i = 1 : size(wrapLocations, 1)
%     item  = (R_transfer* wrapLocations(i,:)')';
%     wrapLocationsRot = [wrapLocationsRot; item];
% end
% 
% SEL_epiShaftRot = (R_transfer * SEL_epiShaft')';
% SEL_pointShaftRot = (R_transfer * SEL_pointShaft')';
% headShaftRot = (R_transfer*headShaft')';
% isthmusShaftRot = (R_transfer*isthmusShaft')';
% % the muscle attachments
% femurShaftLocRotMA = [];
% for i = 1:size(femurMuscle_start,1)
%     item  = (R_transfer* femurShaftLocMA(i,:)')';
%     femurShaftLocRotMA = [femurShaftLocRotMA; item];
% end
% 
% % the markers
% femurShaftLocRotMarkers = [];
% for i = 1:size(markerFemur_start,1)
%     item  = (R_transfer* femurShaftLocMarkers(i,:)')';
%     femurShaftLocRotMarkers = [femurShaftLocRotMarkers; item];
% end
% 
% % figure
% % % scatter3(femurShaftLoc(:,1),femurShaftLoc(:,2),femurShaftLoc(:,3),'black'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % scatter3(femurShaftLocRot(:,1),femurShaftLocRot(:,2),femurShaftLocRot(:,3),'filled'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % % scatter3(femurShaftLocRotMA(:,1),femurShaftLocRotMA(:,2),femurShaftLocRotMA(:,3),'blue'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % % scatter3(femurShaftLocRotMarkers(:,1),femurShaftLocRotMarkers(:,2),femurShaftLocRotMarkers(:,3),'blue'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % % scatter3(SEL_epiShaft(:,1),SEL_epiShaft(:,2),SEL_epiShaft(:,3),'green');
% % scatter3(SEL_epiShaftRot(:,1),SEL_epiShaftRot(:,2),SEL_epiShaftRot(:,3),'red');
% % % scatter3(SEL_pointShaft(:,1),SEL_pointShaft(:,2),SEL_pointShaft(:,3),'green');
% % scatter3(SEL_pointShaftRot(:,1),SEL_pointShaftRot(:,2),SEL_pointShaftRot(:,3),'red');
% % % scatter3(headShaft(:,1),headShaft(:,2),headShaft(:,3),'green');
% % scatter3(headShaftRot(:,1),headShaftRot(:,2),headShaftRot(:,3),'blue');
% % scatter3(isthmusShaftRot(1),isthmusShaftRot(2),isthmusShaftRot(3),'blue')
% % plot3([headShaftRot(:,1),isthmusShaftRot(1)],[headShaftRot(:,2),isthmusShaftRot(2)],[headShaftRot(:,3),isthmusShaftRot(3)],'red');
% % % plot3([SEL_epiShaft(:,1),SEL_pointShaft(:,1)],[SEL_epiShaft(:,2),SEL_pointShaft(:,2)],[SEL_epiShaft(:,3),SEL_pointShaft(:,3)], 'red', 'Linewidth',1.5) % z-axis knee joint
% % plot3([SEL_epiShaftRot(:,1),SEL_pointShaftRot(:,1)],[SEL_epiShaftRot(:,2),SEL_pointShaftRot(:,2)],[SEL_epiShaftRot(:,3),SEL_pointShaftRot(:,3)], 'red', 'Linewidth',1.5) % z-axis knee joint
% % % plot3([0,0],[0,0],[-0.4,0.1], 'black--', 'Linewidth',1.5) % z-axis knee joint
% % % plot3([headShaft(:,1),0],[headShaft(:,2),0],[headShaft(:,3),0], 'black', 'Linewidth',1.5) % z-axis knee joint
% % % plot3([headShaftRot(:,1),0],[headShaftRot(:,2),0],[headShaftRot(:,3),0], 'red', 'Linewidth',1.5) % z-axis knee joint
% 
% %% Find the rotation boxes of the femur
% % innerBox = Femoral Head and neck;
% % middleBox = lesser SEL_point and proximal part of the Shaft;
% % outerBox = Femur proximal to the condylar;
% % Condylar - does not rotate
% 
% % the bone
% HeadNeck = []; LesserTroc = []; Shaft= [];  Condylar =[];
% FemurShaftAxis = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxis = norm(FemurShaftAxis);
% 
% for i = 1:size(femurShaftLocRot,1)
%     itemVector = femurShaftLocRot(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxis)/magn_FemurShaftAxis;
%     if item <= 0.12*(magn_FemurShaftAxis/16)
%         HeadNeck= [HeadNeck; femurShaftLocRot(i,:)];
%     elseif item  < 1.45*(magn_FemurShaftAxis/16) && item > 0.12*(magn_FemurShaftAxis/16)
%         LesserTroc  = [LesserTroc;femurShaftLocRot(i,:)];
%     elseif item < 14*(magn_FemurShaftAxis/16) && item > 1*(magn_FemurShaftAxis/16) % limit for the shaft %0.395
%         Shaft = [Shaft; femurShaftLocRot(i,:)];
%     else
%         Condylar = [Condylar; femurShaftLocRot(i,:)];
%     end
% end
% 
% % Divide the shaft into proximal and distalt part
% ShaftProx = []; ShaftDist = [];
% FemurShaftAxis = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxis = norm(FemurShaftAxis);
% for i = 1:size(Shaft,1)
%     itemVector = Shaft(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxis)/magn_FemurShaftAxis;
%     if item <= 0.5*(magn_FemurShaftAxis/2) % kv: changed from 1 to 0.5 -> smoother
%         ShaftProx= [ShaftProx; Shaft(i,:)];
%     else
%         ShaftDist = [ShaftDist; Shaft(i,:)];
%     end
% end
% 
% % Muscle attachments
% HeadNeckMA = []; LesserTrocMA = []; ShaftMA = [];  CondylarMA =[];
% FemurShaftAxisMA = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxisMA = norm(FemurShaftAxisMA);
% for i = 1:size(femurShaftLocRotMA,1)
%     itemVector = femurShaftLocRotMA(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxisMA)/magn_FemurShaftAxisMA;
%     if item <= 0.12*(magn_FemurShaftAxisMA/16)
%         HeadNeckMA= [HeadNeckMA; femurShaftLocRotMA(i,:)];
%     elseif item  < 1.45*(magn_FemurShaftAxisMA/16) && item > 0.12*(magn_FemurShaftAxisMA/16)
%         LesserTrocMA  = [LesserTrocMA;femurShaftLocRotMA(i,:)];
%     elseif item < 14*(magn_FemurShaftAxisMA/16) && item > (magn_FemurShaftAxisMA/16) % limit for the shaft %0.395
%         ShaftMA = [ShaftMA; femurShaftLocRotMA(i,:)];
%     else
%         CondylarMA = [CondylarMA; femurShaftLocRotMA(i,:)];
%     end
% end
% % Divide the shaft into proximal and distalt part
% ShaftProxMA = []; ShaftDistMA = [];
% FemurShaftAxisMA = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxisMA = norm(FemurShaftAxisMA);
% for i = 1:size(ShaftMA,1)
%     itemVector = ShaftMA(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxisMA)/magn_FemurShaftAxisMA;
%     if item <= 0.5*(magn_FemurShaftAxis/2)
%         ShaftProxMA= [ShaftProxMA; ShaftMA(i,:)];
%     else
%         ShaftDistMA = [ShaftDistMA; ShaftMA(i,:)];
%     end
% end
% 
% % The markers
% HeadNeckMarkers = []; LesserTrocMarkers = []; ShaftMarkers= [];  CondylarMarkers =[];
% FemurShaftAxisMarkers = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxisMarkers = norm(FemurShaftAxisMarkers);
% for i = 1:size(femurShaftLocRotMarkers,1)
%     itemVector = femurShaftLocRotMarkers(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxisMarkers)/magn_FemurShaftAxisMarkers;
%     if item <= 0.12*(magn_FemurShaftAxisMarkers/16)
%         HeadNeckMarkers= [HeadNeckMarkers; femurShaftLocRotMarkers(i,:)];
%     elseif item  < 1.45*(magn_FemurShaftAxisMarkers/16) && item > 0.12*(magn_FemurShaftAxisMarkers/16)
%         LesserTrocMarkers  = [LesserTrocMarkers;femurShaftLocRotMarkers(i,:)];
%     elseif item < 14*(magn_FemurShaftAxisMarkers/16) && item > (magn_FemurShaftAxisMarkers/16) % limit for the shaft %0.395
%         ShaftMarkers = [ShaftMarkers; femurShaftLocRotMarkers(i,:)];
%     else
%         CondylarMarkers = [CondylarMarkers; femurShaftLocRotMarkers(i,:)];
%     end
% end
% % Divide the shaft into proximal and distalt part
% ShaftProxMarkers = []; ShaftDistMarkers = [];
% FemurShaftAxisMarkers = SEL_epiShaftRot - SEL_pointShaftRot; %vector from bottom to top
% magn_FemurShaftAxisMarkers = norm(FemurShaftAxisMarkers);
% for i = 1:size(ShaftMarkers,1)
%     itemVector = ShaftMarkers(i,:)-SEL_pointShaftRot; % vector from each point to the max point
%     %the projection of vector each vector on the largest vector
%     item = dot(itemVector,FemurShaftAxisMarkers)/magn_FemurShaftAxisMarkers;
%     if item <= 0.5*(magn_FemurShaftAxisMarkers/2)
%         ShaftProxMarkers= [ShaftProxMarkers; ShaftMarkers(i,:)];
%     else
%         ShaftDistMarkers = [ShaftDistMarkers; ShaftMarkers(i,:)];
%     end
% end
% 
% % figure
% % scatter3(HeadNeck(:,1),HeadNeck(:,2),HeadNeck(:,3),'black'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % % scatter3(HeadNeckMA(:,1),HeadNeckMA(:,2),HeadNeckMA(:,3),'red');
% % % scatter3(HeadNeckMarkers(:,1),HeadNeckMarkers(:,2),HeadNeckMarkers(:,3),'red');
% % scatter3(LesserTroc(:,1),LesserTroc(:,2),LesserTroc(:,3),'blue');
% % % scatter3(LesserTrocMA(:,1),LesserTrocMA(:,2),LesserTrocMA(:,3),'red');
% % % scatter3(LesserTrocMarkers(:,1),LesserTrocMarkers(:,2),LesserTrocMarkers(:,3),'red');
% % scatter3(Shaft(:,1),Shaft(:,2),Shaft(:,3),'black')
% % scatter3(ShaftMA(:,1),ShaftMA(:,2),ShaftMA(:,3),'red')
% % % scatter3(ShaftMarkers(:,1),ShaftMarkers(:,2),ShaftMarkers(:,3),'red')
% % scatter3(Condylar(:,1),Condylar(:,2),Condylar(:,3),'green')
% % % scatter3(CondylarMA(:,1),CondylarMA(:,2),CondylarMA(:,3),'red')
% % % scatter3(CondylarMarkers(:,1),CondylarMarkers(:,2),CondylarMarkers(:,3),'red')
% % % scatter3(SEL_pointShaftRot(:,1),SEL_pointShaftRot(:,2),SEL_pointShaftRot(:,3),'blue');
% 
% 
% % figure
% % scatter3(ShaftProx(:,1),ShaftProx(:,2),ShaftProx(:,3),'black'); hold on; axis equal; xlabel('x');ylabel('y');zlabel('z')
% % % scatter3(ShaftProxMA(:,1),ShaftProxMA(:,2),ShaftProxMA(:,3),'blue');
% % % scatter3(ShaftProxMarkers(:,1),ShaftProxMarkers(:,2),ShaftProxMarkers(:,3),'blue');
% % scatter3(ShaftDist(:,1),ShaftDist(:,2),ShaftDist(:,3),'red');
% % % scatter3(ShaftDistMA(:,1),ShaftDistMA(:,2),ShaftDistMA(:,3),'blue');
% % % scatter3(ShaftDistMarkers(:,1),ShaftDistMarkers(:,2),ShaftDistMarkers(:,3),'blue');
% 
% %
% innerBox = [HeadNeck];
% innerBoxMA = [HeadNeckMA]; % the top part of the femur (femoral head and greater SEL_point)
% innerBoxMarker = [HeadNeckMarkers];
% middleBox = [LesserTroc;ShaftProx]; % the lesser SEL_point and proximal part of the femur
% middleBoxMA = [LesserTrocMA; ShaftProxMA];
% middleBoxMarker = [LesserTrocMA; ShaftProxMA];
% outerBox = [HeadNeck; LesserTroc; Shaft]; % everything except for the condylar
% outerBox_less = [ShaftDist;Condylar]; % everything apart form inner and middle box
% 
% end
% 
% function [markerCalcn, markerTibia, markerCalcnNR, markerTibiaNR, markerFemur, markerFemurNR] = OpenSimMarkers(markerset, answerLeg, rightbone)
% 
% %find the markers for each bone for the lower extremities
% markerSize = size(markerset.OpenSimDocument.MarkerSet.objects.Marker,2);
% markerCalcn = []; markerTibia = []; markerCalcnNR = []; markerTibiaNR = []; markerFemur = []; markerFemurNR = [];
% if strcmp(answerLeg, rightbone) == 1;
%     for i = 1:markerSize
%         try
%             bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.body.Text;        
%         catch ME
%             if contains(ME.message, 'body')
%                 bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1, i}.socket_parent_frame.Text;
%             else
%                 error('Unknown Error in OpenSimMarkers regarding bonepart, check try-catch blocks')
%             end
%         end
%         
%         if contains(bonepart, 'calcn_r') == 1
%             markerCalcn = [markerCalcn; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerCalcnNR = [markerCalcnNR; i];
%         elseif contains(bonepart, 'tibia_r') == 1
%             markerTibia = [markerTibia; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerTibiaNR = [markerTibiaNR; i];
%         elseif contains(bonepart, 'femur_r') == 1
%             markerFemur = [markerFemur; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerFemurNR = [markerFemurNR; i];
%         end
%     end
% else
%     for i = 1:markerSize
%         try
%             bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.body.Text;        
%         catch ME
%             if contains(ME.message,'body')
%                 bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1, i}.socket_parent_frame.Text;
%             else
%                 error('Unknown Error in OpenSimMarkers regarding bonepart, check try-catch blocks')
%             end
%         end
%         
%         if contains(bonepart, 'calcn_l') == 1
%             markerCalcn = [markerCalcn; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerCalcnNR = [markerCalcnNR; i];
%         elseif contains(bonepart, 'tibia_l') == 1
%             markerTibia = [markerTibia; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerTibiaNR = [markerTibiaNR; i];
%         elseif contains(bonepart, 'femur_l') == 1
%             markerFemur = [markerFemur; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
%             markerFemurNR = [markerFemurNR; i];
%         end
%     end
% end
% end
% 
% function [output]=coordinatesCorrection(input)
% if isempty(input)
%     output = input;
% else
%     Rx = [1 0 0; 0 cos(pi/2) -sin(pi/2); 0 sin(pi/2) cos(pi/2)];
%     Rzz = [cos(-pi/2) -sin(-pi/2) 0; sin(-pi/2) cos(-pi/2) 0; 0 0 1];
% 
%     R =Rzz*Rx;
%     output = (R*input')';
% end
% end

