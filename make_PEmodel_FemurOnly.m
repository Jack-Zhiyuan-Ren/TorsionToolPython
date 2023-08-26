function [ ready ] = make_PEmodel_FemurOnly( answerModel, deformed_model, answerMarkerSet, deform_bone, which_leg, angle, angle_NS, geometryFolder)
place = [cd '\DEFORMED_MODEL\'];

% what model you want to deform
if strcmp(which_leg, 'R') == 1 &&  strcmp(deform_bone, 'F') == 1;
        answerModel_tmp = [ answerModel];
        answerMarkerSet_tmp = [ answerMarkerSet];
else
    answerModel_tmp = [place answerModel];
    answerMarkerSet_tmp = [place answerMarkerSet];
end
dataModel = xml2struct(answerModel_tmp);

if ~exist("DEFORMED_MODEL\Geometry", 'dir')
    mkdir("DEFORMED_MODEL\Geometry")
end

if strcmp(deform_bone, 'F') && strcmp(which_leg, 'R')
    for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry, 'Mesh')
                if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                        .Mesh, 2) > 1
                    for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                            .Mesh, 2)
                        try
                            vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                                .Mesh{j}.mesh_file.Text;
                            copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                        end
                    end
                else
                    try
                        vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                            .Mesh.mesh_file.Text;
                        copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                    end
                end
            end
        else
            try
                if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                        .GeometrySet.objects.DisplayGeometry, 2) > 1
                    for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                        .GeometrySet.objects.DisplayGeometry, 2)
                        try
                            vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                                .GeometrySet.objects.DisplayGeometry{j}.geometry_file.Text;
                            copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                        end
                    end

                else
                    try
                        vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                            .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
                        copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                    end
                end
            end
        end

        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'components')
            if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components, 'PhysicalOffsetFrame')
                if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame, 'attached_geometry')
                    if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry, 'Mesh')
                        if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                .Mesh, 2) > 1
                            for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                    .Mesh, 2)
                                try
                                    vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                        .Mesh{j}.mesh_file.Text;
                                    copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                                end
                            end
                        else
                            try
                                vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                    .Mesh.mesh_file.Text;
                                copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                            end
                        end
                    end
                end
            end
        end
    end
end
% what you want to name the deformed model
answerNameModel = deformed_model;

% the marker set for this model.
markerset = xml2struct(answerMarkerSet_tmp);
answerLegFemur = which_leg;
answerDegFemur = angle;

for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['femur_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh.mesh_file.Text;
        else
            femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
        end
    end
end

if strcmp(answerLegFemur, 'R') == 1; % Rotation of the right foot
%         FA_preAngle = 17.6;
%         NS_preAngle = 123.3;
        % The added anteversion angle is definded
        angleCorrection = answerDegFemur;% - FA_preAngle;
        FA_angle = -(angleCorrection*(pi/180));
        NS_angle = -((angle_NS)*(pi/180));
        % The geomerty of the right femur is imported
        dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
    else % Rotation of the left foot
%         FA_preAngle = 17.6;
%         NS_preAngle = 123.3;
        % The added anteversion angle is definded
        angleCorrection = answerDegFemur;% - FA_preAngle;
        FA_angle = angleCorrection*(pi/180);
        NS_angle = ((angle_NS)*(pi/180));
        % The geometry of the left femur is imported
        dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
end
    % the script for the rotation of the femur is called.
    femur_ns(dataModel, markerset, answerLegFemur, 'R', FA_angle, NS_angle,...
        answerNameModel,answerMarkerSet, dataFemur, place);

for i = 1 : numel(GeometryPathStarts)
    section = filetext(GeometryPathStarts(i) : GeometryPathEnds(i));
    i1 = strfind(section, '<objects');
    i2 = strfind(section, '</objects');
    if ~isempty(i2)
        section = section(i1:i2);
        nameIdx = strfind(section, 'name=');
        absStart = 0;
        absEnd = 0;
        names = [];
        for j = 1 : numel(nameIdx)
            sep = section(nameIdx(j) + 5);
            sep2 = strfind(section(nameIdx(j) + 6 : end), sep);
            sep2 = sep2(1);
            names{j} = section(nameIdx(j) + 6: nameIdx(j) + 4 + sep2);
            tmp1 = strfind(section(nameIdx(j) - 25 : end), '<');
            tmp1 = nameIdx(j) - 25 + tmp1(1);
            objType = section(tmp1 : nameIdx(j)-2);
            endIdx = strfind(section(tmp1 : end), ['</' objType]);
            endIdx = tmp1 + endIdx;
            section(tmp1-1 : endIdx + length(objType) +1);
            fullText{j} = section(tmp1-1 : endIdx + length(objType) +1);
            if j == 1
                absStart = tmp1-1;
            end
            if j == numel(nameIdx)
                absEnd = endIdx + length(objType) +1;
            end
        end
        [~, order] = sort(names);
        sectionNewOrder = '';
        for j = 1 : numel(order)
            sectionNewOrder = [sectionNewOrder fullText{order(j)}];
        end
        newFileText = strrep(newFileText, section(absStart : absEnd), sectionNewOrder);
    end
end

fid = fopen(['DEFORMED_MODEL/' answerNameModel '.osim'],'w');
fprintf(fid, newFileText);
fclose(fid);
