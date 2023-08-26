function [femurMuscle, femurPlace1, femurNR, femurMuscleType] = femur_MA(dataModel, answerLeg, rightbone)

%% Find the muscles in the model
muscleType = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);
muscleType = muscleType{1};
muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);

%Identify the left and right leg
if strcmp(answerLeg, rightbone) == 1;
    femurMA = 'femur_r';
else
    femurMA = 'femur_l';
end

femurMuscle =[];
femurPlace1 = {};
femurNR = [];
femurMuscleType = {};

muscleTypes = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);

for k = 1 : numel(muscleTypes)
    muscleType = muscleTypes{k};
    muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);

    for i = 1:size(muscles,2)
        if isfield(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}, 'GeometryPath')
            AttachmentSize = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}.GeometryPath.PathPointSet.objects);
        else
            AttachmentSize = [];
        end
        
        for j = 1 : numel(AttachmentSize)
            MuscleAttachments = muscles{1,i}.GeometryPath.PathPointSet.objects.(AttachmentSize{j});
    
            if size(MuscleAttachments, 2) == 1 % only one PathPoint of this kind in this muscle (e.g. rect_fem PathPoint in gait2392)
                try
                    CompareStrings1_femur = strcmp(femurMA, MuscleAttachments.body.Text);
                catch ME
                    if contains(ME.message, 'body')
                        CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments.socket_parent_frame.Text,'/bodyset/',''));
                    else
                        error('Unknown Error in femur_MA check try-catch blocks')
                    end
                end
                if CompareStrings1_femur == 1;
                    try
                        femurMuscle = [femurMuscle; str2num(MuscleAttachments.location.Text)];
                    catch
                        femurMuscle = [femurMuscle; str2num(MuscleAttachments.socket_parent_frame.Text)];
                    end
                    femurMuscleType = [femurMuscleType; muscleType];
                    femurNR = [femurNR; i];
                    place1 = AttachmentSize{j};
                    femurPlace1 = [femurPlace1;place1];
                end
            else    % more of a kind - we need to add {1,%d} after the type
                for ii = 1:size(MuscleAttachments,2)
                    try
                        CompareStrings1_femur = strcmp(femurMA, MuscleAttachments{1,ii}.body.Text);
                    catch ME
                        if contains(ME.message, 'body')
                            CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments{1,ii}.socket_parent_frame.Text,'/bodyset/',''));
                        else
                            error('Unknown Error in femur_MA check try-catch blocks')
                        end
                    end
                    if CompareStrings1_femur == 1;
                        try
                            femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.location.Text)];
                        catch
                            femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.socket_parent_frame.Text)];
                        end
                        femurMuscleType = [femurMuscleType; muscleType];
                        femurNR = [femurNR; i];
                        place1 = [AttachmentSize{j} '{1,' num2str(ii) '}'];
                        femurPlace1 = [femurPlace1;place1];
                    end
                end
            end
        end



    end
end