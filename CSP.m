function [W,lambda,A] = CSP(X1_bandpass,X2_bandpass,num_samples)
%   Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    X1_bandpass
    X2_bandpass
    num_samples
end

arguments (Output)
    W
    lambda
    A
end

    fprintf('\n=== COMMON SPATIAL PATTERNS ===\n\n');
    I_mat = eye(num_samples);
    one_T = ones(num_samples,1);

    X1 = (1/sqrt(num_samples)) * X1_bandpass * (I_mat - (one_T * (one_T')));

    X2 = (1/sqrt(num_samples)) * X2_bandpass * (I_mat - (one_T * (one_T')));

    % Compute the covariance matrix of each class
    S1 = cov(X1');   % S1~[C x C]
    S2 = cov(X2');   % S2~[C x C]

    % Solve the eigenvalue problem S1·W = l·S2·W
    [W,L] = eig(S1, S1 + S2);   % Mixing matrix W (spatial filters are columns)
    lambda = diag(L);           % Eigenvalues
    A = (inv(W))';              % Demixing matrix
    fprintf('Dimensions of X1: %d x %d \n', size(X1));
    fprintf('Dimensions of X2: %d x %d \n', size(X2));
    fprintf('Dimensions of W: %d x %d \n', size(W));
    fprintf('Dimensions of λ: %d x %d \n', size(lambda));
    fprintf('Dimensions of A: %d x %d \n', size(A));
    
end
