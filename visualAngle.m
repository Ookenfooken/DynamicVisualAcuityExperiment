function [ degreePixelsV, degreePixelsH ] = visualAngle( sWidth, sHeight, distance )
% Input screen width height and distance from eye
% Output the number of pixels equal to one degree of visual angle,
% vertically and horizontally
%   .

set(0,'units','pixels')  
%Obtains this pixel information
Pix_SS = get(0,'screensize');

pixelsH=Pix_SS(3); %how many pixels are in the width of the screen
pixels_per_meterH=pixelsH./sWidth; %how many pixels per meter

pixelsV=Pix_SS(4); %how many pixels are in the height of the screen
pixels_per_meterV=pixelsV./sHeight; %how many pixels per meter

degreeMeters=tand(1)*distance; %what is one degree of visual angle in meters on the screen
degreePixelsH=degreeMeters*pixels_per_meterH;
degreePixelsV=degreeMeters*pixels_per_meterV;

end

