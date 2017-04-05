function [FILTEREDX,SMOOTHEDX]=kalman_filterpatients(DATA, A, C, Q, R, INITX, INITV)
%[FILTEREDX,SMOOTHEDX]=kalman_filterpatients(DATA_OBS,DATA_RAW, A, C, Q, R, INITX, INITV,sc)
% Goal:
%   run the kalma filter and the kalman smoother 
% Inputs:
%   DATA_OBS: output of kalman_sort_data, n*1 cell array, n = number of patients;
%       DATA_OBS{:,1} stores a num_measurements  by t array: t = number of vists of a patient,
%           the array stores the meausrements longtitudinally.    
%   DATA_RAW: output of kalman_sort_data, n*3 cell array, n = number of patients;
%       DATA_RAW{:,1} stores a num_measurements  by t array: t = number of vists of a patient,
%           the array stores the meausrements longtitudinally.
%       DATA_RAW{:,2} stores patient ids
%       DATA_RAW{:,3} stores visit numbers in a 1 by t matrix
%   A, C, Q, R, INITX, INITV: KF parameters, output of EM algorithm
%       (learn_kalman)
%   sc: scaling factor.
% Output:
%   FILTEREDX:
%       FILTEREDX{:,1} stores a num_measurements  by t array: t = number of vists of a patient,
%           the array stores the filtered values longtitudinally.
%       FILTEREDX{:,2} stores patient ids
%       FILTEREDX{:,3} stores visit numbers in a 1 by t matrix   
%   SMOOTHEDX:
%       SMOOTHEDX{:,1} stores a num_measurements  by t array: t = number of vists of a patient,
%           the array stores the smoothed values longtitudinally.
%       SMOOTHEDX{:,2} stores patient ids
%       SMOOTHEDX{:,3} stores visit numbers in a 1 by t matrix   
% Xiang Liu, 7/19/2012, liuxiang@umich.edu
num_of_patients=length(DATA);
for i = 1 :num_of_patients
    y=DATA{i,3};
    [xfiltered, Vsmooth, VVsmooth, loglik] = kalman_filter(y, A, C, Q, R, INITX, INITV);
    FILTEREDX{i,1}=xfiltered;
    FILTEREDX{i,2}=DATA{i,1};
    FILTEREDX{i,3}=DATA{i,2};
    
end