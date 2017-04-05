%function directly calculates the probability of the logistic regression
%taking value greater than threshold; 

%will need to add a column and row of 0's to the sigma from the KF to add
%the intercept term of the logistic regression

%assume beta includes intercept term
%assume beta is a row vector, mu is a column vector

function P = normP(givenmu, givensigma, beta, threshold)

% taking sigma and mu as given input, adjust for the intercept term
mu=[1;givenmu]; %add intercept term to given mu

[q,r]=size(givensigma);
sigma=zeros(q+1,r+1); %create zero matrix to be filled with cov inputs
sigma(2:q+1,2:r+1)=givensigma; %fill in given sigma

% calculate E(Z)
muZ = beta*mu;

% calculate Var(Z)
sigma2Z = 0;
for i = 1:length(sigma)
   for j = 1:length(sigma)
       sigma2Z = sigma2Z + beta(i)*beta(j)*sigma(i,j);
   end
end

%standardize random variable and calculate P using threshold
P = 1-normcdf((-log(1/threshold-1)-muZ)/sqrt(sigma2Z));