function DVA_pursuit_constant_stimuli
%
% ___________________________________________________________________
%
% Estimates acuity thresholds using a constant but randomized procedure for stimuli
% moving at various velocities and presentation durations
% 
% Cole Shing, 2016, adapted from Dimitrios DVA_pursuit code
% ___________________________________________________________________

% HISTORY
% Dimitrios Palidis, 2015, adapted from SR-Research example code

%% SET Paths
subDir=pwd; %we will create an output folder here, default pwd which is same directory

%% First calculate the number of pixels per degree of angle
sWidth= .427; %input the width of the screen in m
sHeight= .345;%input height of the screen in m
distance= .80; %input distance from screen in m
[~, degreePixelsH]=visualAngle(sWidth, sHeight, distance); % this is pixels/degree in both directions

%% Set parameters of stimuli
XedgeBorder=0.5; %in degrees of visual angle
%reversals=15; %limit of reversals
rampDur=.1;
conditions1=[60; 75; 90]; %STAT, speed of stimuli
ballSize=10;
gaps=[-ballSize,-ballSize,0,0; -ballSize,ballSize,0,0; 0,0,ballSize,ballSize; 0,0,ballSize,-ballSize;];

%General Settings
fixationTime=.5; % minimum duration of fixation cross before stimulus
fixJit=.5; %max of additional random fixation time 
baseblockSize=90; %trials per baseline blocks
numofblocksforBase=1;%amount of blocks baseline
blockSize=150; % trials per block
numofblocks=2; %amount of blocks for trials
%% do conversions to pixels
XedgePix=round(XedgeBorder*degreePixelsH);
conditions=[conditions1 NaN(length(conditions1),1)]; % make a matrix of conditions 1 and NaN
conditions(:,1)=round(conditions(:,1)*degreePixelsH); % times conditions 1 by pixels/degree in the horizontal
%% calculate durations 
Pix_SS = get(0,'screensize'); %get pixels, Pix_SS(3)= horizontal pixels size, Pix_SS(4) = vertical pixels size
%[maxV,mi]=max(conditions(:,1)); %find the fastest one, maxV = max velocity
for i=1:size(conditions(:,1))
    conditions(i,2)=(Pix_SS(3)-2*XedgePix)/conditions(i,1); % duration= screenwidth/speed (all durations equal length of time it takes longest one to traverse screen
end
%% initialize stuff
gaptrials= blockSize*numofblocks; %has to be divisible by 4, used in the actual trials
basetrials=baseblockSize*numofblocksforBase; %used for base trials
ntrials=gaptrials+basetrials; %total trials, has to be even.
directions=[zeros(1,ntrials/2) ones(1,ntrials/2)]; %creates an array of 0s and 1s
directions= directions(randperm(size(directions,2))); %randomizes the arraway of 0s and 1s
gapID=[ones(1,gaptrials/4) ones(1,gaptrials/4)+1 ones(1,gaptrials/4)+2 ones(1,gaptrials/4)+3]; %create 4 gaps ID, then randomize the order
gapID= gapID(randperm(size(gapID,2)));
gapsize=[ones(1,gaptrials/5) ones(1,gaptrials/5)+1 ones(1,gaptrials/5)+2 ones(1,gaptrials/5)+3 ones(1,gaptrials/5)+4];
gapsize= gapsize(randperm(size(gapsize,2))); %create gapsize of 1-5 then randomize order, have to be divisible by 5
basespeed=[ones(1,basetrials/3) ones(1,basetrials/3)+1 ones(1,basetrials/3)+2]; %speed conditions based on 3 conditions, equal amounts
basespeed= basespeed(randperm(size(basespeed,2)));                              %of each condition for both base trials and normal trials
trialCond=[ones(1,gaptrials/3) ones(1,gaptrials/3)+1 ones(1,gaptrials/3)+2];
trialCond= trialCond(randperm(size(trialCond,2)));
%Returns true if the script is running under GNU/Octave
if ~IsOctave 
    commandwindow;
else
    more off;
end

dummymode = 0;

try
    %%%%%%%%%%
    % STEP 1 %
    %%%%%%%%%%
    
    % Added a dialog box to set your own EDF file name before opening
    % experiment graphics. Make sure the entered EDF file name is 1 to 8
    % characters in length and only numbers or letters are allowed.
        
    if IsOctave %Returns true if the script is running under GNU/Octave
        Outputname = 'DEMO';
    else    
        prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
        dlg_title = 'Create EDF file';
        num_lines= 1;
        def     = {'DEMO'}; %default is DEMO
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        Outputname = strcat(answer{1},'_',datestr(now,'yyyymmdd_HHMMSS'));
        fprintf('Outputname: %s\n', Outputname );
        mkdir(subDir,Outputname);
        outputDir=strcat(subDir,'\',Outputname);
    end
    
    %open txt file to log stimulus paramteters.
    fileID = fopen(strcat(Outputname,'\','SLog_',Outputname),'w');
    fprintf(fileID,'%4s %4s %8s %7s %3s %3s %3s %3s\n','trial','velocity','duration','gapsize','direction','gap','responseR','correct');
    
    %%%%%%%%%%
    % STEP 2 %
    %%%%%%%%%%
    
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.
    screenNumber=max(Screen('Screens'));
    [window, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2); %#ok<*NASGU>
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [winWidth, winHeight] = WindowSize(window);
    Screen('TextSize', window,100);
    
    % define center
    center=[winWidth/2 winHeight/2];
 
    %%%%%%%%%%
    % STEP 3 %
    %%%%%%%%%%
    
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    
    el=EyelinkInitDefaults(window);
    tool = tool_functions(el);
    % We are changing calibration to match task background and target
    % this eliminates affects of changes in luminosity between screens
    % no sound and smaller targets
    el.targetbeep = 0;
    el.backgroundcolour = WhiteIndex(el.window);
    el.calibrationtargetcolour= [0 0 0];
    % for lower resolutions you might have to play around with these values
    % a little. If you would like to draw larger targets on lower res
    % settings please edit PsychEyelinkDispatchCallback.m and see comments
    % in the EyelinkDrawCalibrationTarget function
    el.calibrationtargetsize= 2;
    el.calibrationtargetwidth=0.5;
    % call this function for changes to the el calibration structure to take
    % affect
    EyelinkUpdateDefaults(el);
    
    %%%%%%%%%%
    % STEP 4 %
    %%%%%%%%%%
    
    % Initialization of the connection with the Eyelink tracker
    % exit program if this fails.
    
    if ~EyelinkInit(dummymode)
        %fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
%     % open file to record data to
%     res = Eyelink('Openfile', edfFile);
%     if res~=0
%         fprintf('Cannot create EDF file ''%s'' ', edffilename);
%         cleanup;
%         return;
%     end
%     
%     % make sure we're still connected.
%     if Eyelink('IsConnected')~=1 && ~dummymode
%         cleanup;
%         return;
%     end
    
    %%%%%%%%%%
    % STEP 5 %
    %%%%%%%%%%
    
    % SET UP TRACKER CONFIGURATION
    
    Eyelink('command', 'add_file_preamble_text ''Recorded by DVA1.0''');
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    
    % This command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');
    Eyelink('command', 'generate_default_targets = YES');
    
    % STEP 5.1 retrieve tracker version and tracker software version
    [v,vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    
    if v == 3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        
       % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
    % allow to use the big button on the eyelink gamepad to accept the
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
       
    %%%%%%%%%%
    % STEP 6 %
    %%%%%%%%%%
    
    % Hide the mouse cursor
    Screen('HideCursorHelper', window);
    % enter Eyetracker camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
    
    %%%%%%%%%%
    % STEP 7 %
    %%%%%%%%%%
    
    % Now starts running individual trials
    % You can keep the rest of the code except for the implementation
    % of graphics and event monitoring
    % Each trial should have a pair of "StartRecording" and "StopRecording"
    % calls as well integration messages to the data file (message to mark
    % the time of critical events and the image/interest area/condition
    % information for the trial)
    
    %% BEGIN TRIALS
    
    % Initialize counters and data storage
    breakIn=baseblockSize; %count down to break
    blockcount=1;
    %ii=ones(1,size(conditions,1)); %condition counter
    correct=zeros(1,size(conditions,1)); %track correct in a row for each condition
    %ascending=zeros(1,size(conditions,1)); % tracks if we are ascending staircase for each condition,1 means ascending 0 descending
    %ascending(:)=3; %start out neither ascending or descending
    for i=1:ntrials %go until we have completes reversals for all conditions
    %for j=1:20
        j = i - basetrials; %this is for trials after base trials
        %%%Open EDF file for recording%%%
        disp(strcat(num2str(i),'begin trial'))
        % open file to record data to
        trnum=sprintf('%3.3d',i);
        trialFile=strcat(trnum); %name of edf with trail index
        res = Eyelink('Openfile',trialFile);
        if res~=0
            %fprintf('Cannot create EDF file ''%s'' ', edffilename);
            cleanup;
            return;
        end
        disp(strcat(num2str(i),'file open'))
        % make sure we're still connected.
        if Eyelink('IsConnected')~=1 && ~dummymode
            cleanup;
            return;
        end    
        
        Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        if breakIn==0; % enter break/calibration
            if blockcount == numofblocksforBase %after the blockcount equals base number of blocks, switch to trials
               breakIn=blockSize; 
            else
               blockcount=blockcount+1; %checking blockcount 
               breakIn=baseblockSize; %renew breakIn to full blocksize for baseline
            end
            EyelinkDoTrackerSetup(el); %enable calibration/drift correction
        end
        
        % choose random trial condition and get parameters
        %[trialCond,speed,duration]=prepareRandomTrial(conditions,revs,reversals);
        if i > basetrials
            speed=conditions(trialCond(1,j),1);
            duration=conditions(trialCond(1,j),2);
            conditionused=conditions1(trialCond(1,j),1);
        else
            speed=conditions(basespeed(1,i),1);
            duration=conditions(basespeed(1,i),2);
            conditionused=conditions1(basespeed(1,i),1);
        end
            
        % write trial parameters to logfile         
        
        % STEP 7.1
        % Sending a 'TRIALID' message to mark the start of a trial in Data
        % Viewer.  This is different than the start of recording message
        % START that is logged when the trial recording begins. The viewer
        %will not parse any messages, events, or samples, that exist in
        % the data file prior to this message.
        Eyelink('Message', 'TRIALID %d', i);
        
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', i,8);
        % Before recording, we place reference graphics on the host display
        % Must be in offline mode to transfer image to Host PC
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen %d', 0);
         
        
        if directions(i)==1 %Set cross starting position
            cross=[XedgePix+speed*rampDur center(2)]; %if going left to right
        else
            cross=[winWidth-XedgePix-speed*rampDur center(2)]; % if going right to left
        end
        
        if directions(i)==1 % set starting position for the ball
            xo =  cross(1)-speed*rampDur; %going left to right         
            yo =  center(2);
        else
            xo =  cross(1)+speed*rampDur; %going right to left
            yo =  center(2);
        end        
        
        x=xo;
        y=yo;
        ball([1 3]) = [x-ballSize x+ballSize]; %create ball size
        ball([2 4]) = [y-ballSize y+ballSize];
        WaitSecs(0.1);
        % STEP 7.2
        % Do a drift correction at the beginning of each trial
        % Performing drift correction (checking) is optional for
        % EyeLink 1000 eye trackers. Drift correcting at different
        % locations x and y depending on where the ball will start
        % we change the location of the drift correction to match that of
        % the target start position
        % Note drift correction does not accept fractionals in PTB!
       % EyelinkDoDriftCorrection(el,round(x),round(y));
        
        % STEP 7.3
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
        
        % get eye that's tracked
        eye_used = Eyelink('EyeAvailable');
        
        
        %%%display fixation cross for fixationtime seconds
        
        %Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('FillRect', window, el.backgroundcolour);
        Screen('DrawLines',window,[-10 10 0 0; 0 0 -10 10],1,[0,0,0],cross);
        Screen('Flip', window);
        pos= [x y];
        WaitSecs(fixationTime); %wait for the specified fixation time
        jit=fixJit*rand;
        WaitSecs(jit);%wait for an additional random fixation time
        sttime = GetSecs; 
        trialTime = sttime + duration; 
        
%         stim_onseted=0;
        
        while GetSecs < trialTime
                
            % STEP 7.4
            % Prepare and show the screen.
            % Enable alpha blending with proper blend-function. We need it
            % for drawing of smoothed points:
            %Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('FillRect', window, el.backgroundcolour);
            %
            
            time = (GetSecs-sttime);
             %if the duration of step ramp elapsed and in trials draw gap, otherwise draw circle, if duration not elapsed draw cross
             if  i > basetrials
                 Screen('FrameOval', window,[0 0 0], ball, 3);
                 if time>= rampDur
                     Screen('DrawLine', window, el.backgroundcolour, x+gaps(gapID(j),1), y+gaps(gapID(j),2), x+gaps(gapID(j),3), y+gaps(gapID(j),4),gapsize(j));
                 end
             elseif i <= basetrials %time>= rampDur && 
                 tool.draw_target_gaussian(ball, [0 0 0]);
             end
            
            Screen('Flip', window);
%             if stim_onseted==0;
%                 Eyelink('Message', 'Stim_onset');
%                 stim_onset=1;
%             end
            Eyelink('Message', 'SYNCTIME');
            % STEP 7.5
            % send the location of the target at each iteration so that
            % target can be displayed in Dataviewer
            Eyelink('message', '!V TARGET_POS TARG1 (%d, %d) 1 0',floor(x),floor(y));
            
            %y = y;
            if (directions(i)==1) %move x along the screen
                x =  xo + speed*time;
            else
                x =  xo - speed*time;
            end
                                
            ball([1 3]) = [x-ballSize x+ballSize]; %update ball location
            ball([2 4]) = [y-ballSize y+ballSize];
        end
      
        % STEP 7.6
        % add 100 msec of data to catch final events and blank display
        Screen('FillRect', window, el.backgroundcolour);
        
        Screen('Flip', window);
        WaitSecs(0.1);
        Eyelink('StopRecording');
        
        % close and save edf file for trial
        Eyelink('Command', 'set_idle_mode');
        
        if i > basetrials
            %acceptedInputs=[36 35 34 33 81 67]; % what keyboaard input codes are accepted (must test command to see what number corresponds to what key)
            acceptedInputs=[103 97 99 105 81 67];%macbook q,a,w,s,z,c
            % these numbers correspond to: 7,1,3,9,q,c (numbers are on number keypad
            p=keyboardInput(acceptedInputs);
            inputMapping=acceptedInputs;
            response=find(inputMapping==p); %stuck here till accepted input is pressed
        end
        WaitSecs(0.2);
        Eyelink('CloseFile'); %done recording file
       
    
        try
            %fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile',trialFile,outputDir,1);
            %status=Eyelink('ReceiveFile',trialFile,strcat(subDir,edfFile),'dest_is_path');
            %movefile(strcat(edfFile,num2str(i)),strcat(pwd,edfFile));
            %if status > 0
               % fprintf('ReceiveFile status %d\n', status);
            %end
            %if 2==exist(edfFile, 'file')
                %fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            %end
        catch %#ok<*CTCH>
            %fprintf('Problem receiving data file ''%s''\n', edfFile );
        end
        disp(strcat(num2str(i),'recieved file'))
        
        if i > basetrials
            if response==gapID(j) %check if they were correct
                trialCorrect=1;
            else
                trialCorrect=0;
            end
            gapinput=gapID(j);
            gapspace=gapsize(j);
        else
            response=0;
            trialCorrect=-1;
            gapinput=0;
            gapspace=0;
        end
        %fprintf(fileID,'%1d %6d %1d %2d %5d',i, conditions1(trialCond,1), duration, gapsize(trialCond),directions(i));       
        fprintf(fileID,'%1d %6d %1d %2d %5d %10d %12d %9d \n',i, conditionused, duration, gapspace,directions(i),gapinput, response, trialCorrect);
        
        Screen('FillRect', window, el.backgroundcolour);
        Screen('Flip', window);
        WaitSecs(0.5); % inter-trial-interval
        
        if response==5 %input 'q' to quit
             close all
             error('escape command issued')
        elseif response==6 %input 'c' to calibrate
            EyelinkDoTrackerSetup(el);
%         elseif (fixationFailed(fixbound,fixationStart, degreePixelsH, degreePixelsV, center,eye_used)==1) %if we failed to maintain fixation in radius, then show message and skip staircase modifications
%             Screen('FillRect', window, [999 999 999]);
%             Screen('DrawText',window,'FIXATE!!!',center(1)-600, center(2)-100,[999 999 999] , [0 0 1]);
%             Screen('Flip', window);
%             WaitSecs(0.7);
%             fprintf(fileID,'%10d %12d %9d\n',1, 3, 3);        
        end
           % ii(trialCond)=ii(trialCond)+1;
            breakIn=breakIn-1; %minus one to breakin for blocksize              
        
        % STEP 7.7
        % Send out necessary integration messages for data analysis
        % See "Protocol for EyeLink Data to Viewer Integration-> Interest
        % Area Commands" section of the EyeLink Data Viewer User Manual
        % IMPORTANT! Don't send too many messages in a very short period of
        % time or the EyeLink tracker may not be able to write them all
        % to the EDF file.
        % Consider adding a short delay every few messages.
        WaitSecs(0.001);
        % Send messages to report trial condition information
        % Each message may be a pair of trial condition variable and its
        % corresponding value follwing the '!V TRIAL_VAR' token message
        % See "Protocol for EyeLink Data to Viewer Integration-> Trial
        % Message Commands" section of the EyeLink Data Viewer User Manual
        WaitSecs(0.001);
        
        
        Eyelink('Message', '!V TRIAL_VAR index %d', i);
        
        % a limitation of the currect ETB only accepts ints as input to
        % messages and commands a possible work around is given below
        
        
%         msg1 = sprintf('!V TRIAL_VAR freq_x %2.3f ', trials(1,i));
%         msg2 = sprintf('!V TRIAL_VAR freq_y %2.3f ', trials(2,i));
%         Eyelink('Message', msg1);
%         Eyelink('Message', msg2);     
        
        % STEP 7.8
        % Sending a 'TRIAL_RESULT' message to mark the end of a trial in
        % Data Viewer. This is different than the end of recording message
        % END that is logged when the trial recording ends. The viewer will
        % not parse any messages, events, or samples that exist in the data
        % file after this message.
        Eyelink('Message', 'TRIAL_RESULT 0');
    i=i+1;  
    %log=[log;[i conditions1(trialCond,1) duration eccentricity gapsize(trialCond) leftRight(i) directions(i)]]; 
    save(strcat(Outputname,'\',Outputname)); %save workspace after each trial
    
    %if escape key is pressed, abort program
    [~, ~, keyCode] = KbCheck;
    if keyCode(KbName('ESCAPE'))        
            KbReleaseWait; 
            disp(strcat(num2str(i),'Escape key pressed. Trial aborted.'))
            break;
    end
    
    end %main for loop end
    
    %%%%%%%%%%
    % STEP 8 %
    %%%%%%%%%%
    
    % End of Experiment; close the file first
    % close graphics window, close data file and shut down tracker
    
   % figure
  %  for i =1:size(conditions1,1)
%        revs1=revPoints(i,2:length(revPoints(i,:)));
%        threshold= NaNmean(revs1);
   %     subplot(4,3,i)
     %   plot(staircase(1,:,i))
     %   hold on
     %   scatter([1:size(staircase(2,:,i),2)],staircase(1,:,i),[],[staircase(2,:,i)' staircase(2,:,i)' staircase(2,:,i)'].*999,'filled','MarkerEdgeColor','k')
     %   hold off
%        title(strcat(edfFile,' StairCase',' ',' V: ',num2str(conditions1(i)),' Threshold: ',' ',num2str(threshold)));
  %  end
    
    
    %%%%%%%%%%
    % STEP 9 %
    %%%%%%%%%%
    
    % run cleanup function (close the eye tracker and window).
        cleanup;
    
catch
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    %Eyelink('CloseFile');
    fclose(fileID); %close stimulus log file
    
%     try
%         %fprintf('Receiving data file ''%s''\n', edfFile );
%         status=Eyelink('ReceiveFile'); %why receieve aother file?
%         if status > 0
%             %fprintf('ReceiveFile status %d\n', status);
%         end
%         if 2==exist(edfFile, 'file')
%             %fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
%         end
%     catch %#ok<*CTCH>
%         %fprintf('Problem receiving data file ''%s''\n', edfFile );
%     end 
    cleanup;
    %fprintf('%s: some error occured\n', mfilename);
    psychrethrow(lasterror); %#ok<*LERR>
    
end

    function cleanup
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        Screen('CloseAll');
    end

end
