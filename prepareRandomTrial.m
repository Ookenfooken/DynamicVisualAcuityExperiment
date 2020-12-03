function [ trialCond,speed,duration ] = prepareRandomTrial( conditions,revs,reversals )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
        trialCond=ceil(rand*size(conditions,1)); %generate from 1 to size of conditions
        while(revs(trialCond)==reversals+1) %does nothing????
            trialCond=ceil(rand*size(conditions,1));
        end
        speed=conditions(trialCond,1);
        duration=conditions(trialCond,2);


end

