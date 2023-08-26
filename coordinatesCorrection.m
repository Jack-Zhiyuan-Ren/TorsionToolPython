function [output]=coordinatesCorrection(input)
if isempty(input)
    output = input;
else
    Rx = [1 0 0; 0 cos(pi/2) -sin(pi/2); 0 sin(pi/2) cos(pi/2)];
    Rzz = [cos(-pi/2) -sin(-pi/2) 0; sin(-pi/2) cos(-pi/2) 0; 0 0 1];

    R =Rzz*Rx;
    output = (R*input')';
end