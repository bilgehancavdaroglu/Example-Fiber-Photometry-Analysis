function [normDat] = deltaFF3 (dat1,BinSize,Window)
% DELTAFF3
% Function to normalize data.
% Inputs:
%   dat1: Raw channel data
% Output:
%   normDat: Normalized data

%This is closer to what is written in doric manual and uses moving average
%
[short,long]=movavg(dat1,1,Window/BinSize,0);%5/BinSize ensures that the window is Window seconds regardless of bin size
%Find deltaF/F: Subtract control from raw data
normDat = (dat1 - long)./ long; %this gives deltaF/F
normDat = normDat * 100; % get %