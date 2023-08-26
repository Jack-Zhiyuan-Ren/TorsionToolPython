function [markerCalcn, markerTibia, markerCalcnNR, markerTibiaNR, markerFemur, markerFemurNR] = OpenSimMarkers(markerset, answerLeg, rightbone)

%find the markers for each bone for the lower extremities
markerSize = size(markerset.OpenSimDocument.MarkerSet.objects.Marker,2);
markerCalcn = []; markerTibia = []; markerCalcnNR = []; markerTibiaNR = []; markerFemur = []; markerFemurNR = [];
if strcmp(answerLeg, rightbone) == 1;
    for i = 1:markerSize
        try
            bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.body.Text;        
        catch ME
            if contains(ME.message, 'body')
                bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1, i}.socket_parent_frame.Text;
            else
                error('Unknown Error in OpenSimMarkers regarding bonepart, check try-catch blocks')
            end
        end
        
        if contains(bonepart, 'calcn_r') == 1
            markerCalcn = [markerCalcn; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerCalcnNR = [markerCalcnNR; i];
        elseif contains(bonepart, 'tibia_r') == 1
            markerTibia = [markerTibia; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerTibiaNR = [markerTibiaNR; i];
        elseif contains(bonepart, 'femur_r') == 1
            markerFemur = [markerFemur; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerFemurNR = [markerFemurNR; i];
        end
    end
else
    for i = 1:markerSize
        try
            bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.body.Text;        
        catch ME
            if contains(ME.message,'body')
                bonepart = markerset.OpenSimDocument.MarkerSet.objects.Marker{1, i}.socket_parent_frame.Text;
            else
                error('Unknown Error in OpenSimMarkers regarding bonepart, check try-catch blocks')
            end
        end
        
        if contains(bonepart, 'calcn_l') == 1
            markerCalcn = [markerCalcn; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerCalcnNR = [markerCalcnNR; i];
        elseif contains(bonepart, 'tibia_l') == 1
            markerTibia = [markerTibia; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerTibiaNR = [markerTibiaNR; i];
        elseif contains(bonepart, 'femur_l') == 1
            markerFemur = [markerFemur; str2num(markerset.OpenSimDocument.MarkerSet.objects.Marker{1,i}.location.Text)];
            markerFemurNR = [markerFemurNR; i];
        end
    end
end
