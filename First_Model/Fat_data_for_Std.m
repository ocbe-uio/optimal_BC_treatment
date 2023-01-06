%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate computations related to fat
%
% Author: Tugba Akman Date: Dec 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are taken from
% Obesity-Activated Adipose-Derived Stromal Cells Promote
% Breast Cancer Growth and Invasion
% Neoplasia (2018) 20, 1161–1174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fat_CD_ydata,Fat_HFD_ydata]=Fat_data_for_Std

tumor_CD_ydata = [138.92, 324.7, 480.76]';
tumor_HFD_ydata = [326.8, 702.025, 1114.833]';

diam_CD = 0.1;
diam_HFD = 0.2;
radius_CD=diam_CD/2;
radius_HFD=diam_HFD/2;

ASC_data_matrix_CD = [14, 36, 35, 41];
ASC_data_matrix_HFD = [27, 48, 50, 32, 18];

w_CD = ASC_data_matrix_CD./(0.22004866);

w_HFD = ASC_data_matrix_HFD./(0.22004866);

for i=1:4
    NumbOfASC_CD_ydata = ((tumor_CD_ydata(3))^(2/3))*(w_CD(i)/diam_CD); %143.18 form the image
    Fat_CD_ydata(i) = NumbOfASC_CD_ydata*(4/3)*pi*(radius_CD^3);
end

for i=1:5
    NumbOfASC_HFD_ydata = ((tumor_HFD_ydata(3))^(2/3))*(w_HFD(i)/diam_HFD); %159.09 from the image
    Fat_HFD_ydata(i) = NumbOfASC_HFD_ydata*(4/3)*pi*(radius_HFD^3); %4.6013e+05
end

end
