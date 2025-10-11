function [X1_CSP,X2_CSP] = CSP(X1_bandpass,X2_bandpass,num_samples)

%  Maximizes the variance of the spatially filtered signal under one condition while minimizing it for the other.

arguments (Input)
    X1_bandpass
    X2_bandpass
    num_samples
end

arguments (Output)
    X1_CSP
    X2_CSP
end

    fprintf('\n=== COMMON SPATIAL PATTERNS ===\n\n');
    I_mat = eye(num_samples);
    one_T = ones(num_samples,1);

    X1_mean_centered = (1/sqrt(num_samples)) * X1_bandpass * (I_mat - (one_T * (one_T')));

    X2_mean_centered = (1/sqrt(num_samples)) * X2_bandpass * (I_mat - (one_T * (one_T')));

    % Compute the covariance matrix of each class
    S1 = cov(X1_mean_centered');
    S2 = cov(X2_mean_centered');

    % Solve the eigenvalue problem S1·W = l·S2·W
    [W,L] = eig(S1, S1 + S2);   % Filter matrix W
    lambda = diag(L);           % Eigenvalues λ
    [lambda, idx] = sort(lambda, 'descend');
    % Sort the eigenvalues and corresponding eigenvectors
    W = W(:, idx);
    A = (inv(W))';              % Pattern matrix A
    X1_CSP = (W(: , 1:2))' * X1_mean_centered;
    X2_CSP = (W(: , 1:2))' * X2_mean_centered;

    fprintf('Dimensions of X1_CSP: %d x %d \n', size(X1_CSP));
    fprintf('Dimensions of X2_CSP: %d x %d \n', size(X2_CSP));
    fprintf('Dimensions of W_CSP: %d x %d \n', size(W));
    fprintf('Dimensions of λ_CSP: %d x %d \n', size(lambda));
    fprintf('Dimensions of A: %d x %d \n', size(A));
    
end
