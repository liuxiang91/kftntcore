function [ A0 C0 Q0 R0 INITX0 INITV0 ] = initializeEM( data )
%INITIALIZEEM Summary of this function goes here
%   Detailed explanation goes here
numPat=size(data,1)-1;
INITX0=zeros(9,1);
for i =1:numPat
    INITX0=INITX0+data{i+1,3}(:,1);
end
INITX0=INITX0/numPat;
[TM, ERR ] = kalman_est_cov( data(2:end,:) );
A0 = TM;
C0 = eye(9,9);
INITV0 = ERR;
Q0 = ERR./2;
R0 = ERR./2;

Q0=triu(tril(Q0));

R0=triu(tril(R0));


end

