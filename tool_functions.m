function tool = tool_functions(el)
    
    %adjustable parameters for gaussian target
    %gaussSize = 40; %99; %% we want this to be 1 deg in size roughly
    gaussSize = 120;
    gaussVar = 8;
    %gaussVar = 80;

%PsychToolbox wrapper for Eyecatch
    tool.draw_target_gaussian        = @(pos, rgb)draw_target_gaussian(pos, rgb);
%     tool.get_center         = @get_center;
%     tool.screen_to_center   = @(pos)screen_to_center(pos);
%     tool.center_to_screen   = @(pos)center_to_screen(pos);    
    
    gaussTexture = zeros(gaussSize, gaussSize, 4); %RGBA
    gaussCenter = [(gaussSize + 1) / 2, (gaussSize + 1) / 2];
    
    for ii = 1:gaussSize
        for jj = 1:gaussSize
            %gaussTexture(ii,jj,4) = -exp(-((ii - gaussCenter(1))^2 + (jj - gaussCenter(1))^2) / (2 * gaussVar^2)) * 127.5 + 127.5;
            gaussTexture(ii,jj,4) = exp(-((ii - gaussCenter(1))^2 + (jj - gaussCenter(1))^2) / (2 * gaussVar^2)) * 255;
        end
    end
    
    gaussIdx = Screen( 'MakeTexture', el.window, gaussTexture);
    
    function draw_target_gaussian(pos, ~)
        %pos_ = center_to_screen(pos); 
        pos_ = pos;
        Screen( 'DrawTexture', el.window, gaussIdx, [], [ pos_(1) - gaussCenter(1), pos_(2) - gaussCenter(2), pos_(1) + gaussCenter(1), pos_(2) + gaussCenter(2) ] );
    end

%     function screenCenter = get_center()
%         [screenWidth, screenHeight] = WindowSize(el.window);
%         screenCenter = [round(screenWidth / 2), round(screenHeight / 2)];
%     end
% 
%     function centerCoord = screen_to_center(screenCoord)
%         centerCoord = screenCoord;
%         centerCoord(1:2) = [screenCoord(1), screenCoord(2)] - get_center();
%         centerCoord(2) = -centerCoord(2);
%     end
% 
%     function screenCoord = center_to_screen(centerCoord)
%         screenCoord = centerCoord;
%         screenCoord(1:2) = [centerCoord(1), -centerCoord(2)] + get_center();
%     end

    
end









