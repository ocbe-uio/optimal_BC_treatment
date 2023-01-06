%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to calculate var and std
%
% Author: Tugba Akman Date: Dec 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stdev1,weight1,stdev2,weight2] = statistical_tools_HFD(data1,data2)

weight1 = var(data1');
weight2 = var(data2');
stdev1 = std(data1');
stdev2 = std(data2');

end