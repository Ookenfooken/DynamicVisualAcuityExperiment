function screenCoord = center_to_screen(centerCoord)
    [winWidth, winHeight] = WindowSize(window);
    % define center
    center = [winWidth/2 winHeight/2];
    screenCoord = centerCoord;
    screenCoord(1:2) = [centerCoord(1), -centerCoord(2)] + center;
end

