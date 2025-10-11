function [X1_LDA, X2_LDA] = LDA(X1_CSP, X2_CSP)

%   Transforms features into a lower dimentional space.

arguments (Input)
    X1_CSP
    X2_CSP
end

arguments (Output)
    X1_LDA
    X2_LDA
end

fprintf('\n=== LINEAR DISCRIMINANT ANALYSIS ===\n\n');

% Perform LDA transformation on the input CSP data
sumX1 = sum(X1_CSP.^2, 2);
sumX2 = sum(X2_CSP.^2, 2);
logX1 = log(sumX1);
logX2 = log(sumX2);

meanX1 = sum(logX1', 1)/2;
meanX2 = sum(logX2', 1)/2;
mean = (meanX1 + meanX2)/2;

Sb = 2*((meanX1 - mean)' * (meanX1 - mean) + (meanX2 - mean)' * (meanX2 - mean)); % Between-class scatter

Sw1 = (logX1' - meanX1)' * (logX1' - meanX1);
Sw2 = (logX2' - meanX2)' * (logX2' - meanX2);
Sw = Sw1 + Sw2; % Within-class scatter

% Compute the LDA transformation matrix
W = Sw\Sb; % W = Inverse (Sw) * Sb

[V,L] = eig(W);
lambda = diag(L);

% Sort the eigenvalues and corresponding eigenvectors
[lambda, idx] = sort(lambda, 'descend');
V = V(:, idx);

% Compute the projection of the CSP data onto the LDA space
X1_LDA = logX1' * V(:, 1); % Project onto the first LDA components
X2_LDA = logX2' * V(:, 1); % Project onto the second LDA components

fprintf('Dimensions of X1_LDA: %d x %d \n', size(X1_LDA));
fprintf('Dimensions of X2_LDA: %d x %d \n', size(X2_LDA));
fprintf('Dimensions of W_LDA: %d x %d \n', size(W));
fprintf('Dimensions of EigenVector matrix: %d x %d \n', size(V));

end
