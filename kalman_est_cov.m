function [TM, COVMAT ] = kalman_est_cov( DATA )
%[TM, COVMAT ] = kalman_est_cov( DATA, num_measures )
% Goals:
%   Calculate the regression coefficient and the covariance of the
%   residuals.
% Inputs:
%   DATA: output of kalman_sort_data, cell array
%       DATA{:,1} stores a num_measurements  by t array: t = number of vists of a patient,
%       the array stores the meausrements longtitudinally.
%   num_measures: the dimension of the states, usually = 9.
% Output:
%   TM: regression coefficient
%   COVMAT:covariance of the residuals
% Xiang Liu, 7/19/2012, liuxiang@umich.edu
numPat = size(DATA,2);


Y = DATA{1,3}(:,2:end);
X = DATA{1,3}(:,1:end-1);

for i = 2:numPat
    Y = [Y DATA{i,3}(:,2:end)];
    X = [X DATA{i,3}(:,1:end-1)];
end

[TM,BINT,ERR] = regress(Y(1,:)',X');
TM = TM';
for i = 2:9
    [B,BINT,R] = regress(Y(i,:)',X');
    TM = [TM; B'];
    ERR = [ERR R];
end
for i = 1:size(ERR,2)
    ERRNEW(:,i)=ERR(~isnan(ERR(:,i)), i);
end
COVMAT = cov(ERRNEW);
end

