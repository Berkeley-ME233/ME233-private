function [Lambda,idx] = me233_autocov(w,y,N)
% [LAMBDA,IDX] = ME233_AUTOCOV(W,Y,N)
%     Calculates auto-covariance for [w y]' for all correlation indices
%     from -N to +N, where the convention used is
%         L(j) = E[x(k+j) x'(k)]
%     The resulting auto-covariance is normalized so that L(0) has 
%     ones along its diagonal elements
%     
% Inputs:
%     W       vector of doubles representing a random ergodic sequence
%     Y       vector of doubles representing a random ergodic sequence
%     N       maximum correlation index to be considered
%     
% Outputs:
%     LAMBDA  ((2N+1) x 4) vector of doubles representing the
%             auto-correlation of x=[w y]'; each row represents a
%             particular value of the correlation index 
%             LAMBDA(:,1) represents E[w(k+j) w(k)] as a function of j
%             LAMBDA(:,2) represents E[y(k+j) w(k)] as a function of j
%             LAMBDA(:,3) represents E[w(k+j) y(k)] as a function of j
%             LAMBDA(:,4) represents E[y(k+j) y(k)] as a function of j            
%     IDX     correlation indices for each row of LAMBDA

%Richard Conway, 2/5/06



numPts = length(w);



%input error checking
if nargin < 3
    error('Insufficient inputs into me233_autocov')
elseif ~isa(w,'double') || ~isvector(w)
    error(['The first input into me233_autocov must be a vector ',...
        'of doubles'])
elseif ~isa(y,'double') || ~isvector(y)
    error(['The second input into me233_autocov must be a vector ',...
        'of doubles'])
elseif ~isa(N,'double') || ~isscalar(N)
    error('The third input into me233_autocov must be a scalar double')
elseif length(y) ~= numPts
    error(['The first and second inputs into me233_autocov must ',...
        'have the same length']);
end
if ~isreal(w) || ~isreal(y) || ~isreal(N)
    warning('Ignoring imaginary parts of inputs')
    x = real(x);
    y = real(y);
    N = real(N);
end



%turns w and y into column vectors
if size(w,2) ~= 1
    w = w';
end
if size(y,2) ~= 1
    y = y';
end

%subtracts mean from w and y and normalizes them 
%with respsect to their covariances
w = w - mean(w);
y = y - mean(y);
w = w / sqrt(cov(w));
y = y / sqrt(cov(y));



%initializes loop variables
Lambda = zeros(2*N+1,4);
j = -N:N;



%calculates auto-covariance and cross-covariance for negative
%correlation coefficient
for i = 1:N
    idx_low = 1:numPts+j(i);
    idx_high = 1-j(i):numPts;
    
    mat = [w(idx_low) w(idx_high) y(idx_low) y(idx_high)];
    temp = cov(mat);
    Lambda(i,1) = temp(1,2);
    Lambda(i,2) = temp(2,3);
    Lambda(i,3) = temp(1,4);
    Lambda(i,4) = temp(3,4);
end



%calculates cross-covariance
temp = cov([w y]);
Lambda(N+1,:) = temp([1 2 3 4]);



%calculates auto-covariance and cross-covariance for positive
%correlation coefficient
for i = N+2:2*N+1
    idx_low = 1:numPts-j(i);
    idx_high = 1+j(i):numPts;
    
    mat = [w(idx_low) w(idx_high) y(idx_low) y(idx_high)];
    temp = cov(mat);
    Lambda(i,1) = temp(1,2);
    Lambda(i,2) = temp(1,4);
    Lambda(i,3) = temp(2,3);
    Lambda(i,4) = temp(3,4);
end

idx = j';
