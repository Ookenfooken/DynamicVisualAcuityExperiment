function [ p ] = keyboardInput( acceptedInputs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
         WaitSecs(.2)
        [s,key,ds]=KbWait;
        p=find(key==1);
        p=p(1);
        %any(find(acceptedInputs==p))
        while ~(any(p==acceptedInputs))
            [secs,key,dsecs]=KbWait;
            p=find(key==1);
            p=p(1);
        end

end

