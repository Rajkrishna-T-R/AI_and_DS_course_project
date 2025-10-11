function [filtered_data_final] = FIR(path,num_trials,high_cutoff,low_cutoff)

% Isolate the brain activity in a specific frequency range
  
arguments (Input)
 path
 num_trials
 high_cutoff
 low_cutoff
end

arguments (Output)
  filtered_data_final
end

fprintf('\n=== FINITE IMPULSE RESPONSE === \n\n')

data=load(path);

all_trials = data.eeg_data_wrt_task_rep_no_eog_256Hz_end_trial;


fs = 256;              % Sampling frequency (Hz)  
order = 100;           % Filter order (higher = sharper) 
num_channels=64;
num_samples=512;


 % More commonly used order is higher order ex:100 for EEG 


% filter design

Wn = [low_cutoff high_cutoff] / (fs/2); % Bandpass low to up Hz, normalized

b = fir1(order, Wn, 'bandpass');        % FIR filter coefficients

filtered_data_final = cell(2,num_trials); % 2 x 100

% each cell is  [64 x 512]

% Filtering Each channel in one trail--

% ----- class 1 -----


for trial_num = 1:num_trials

    trial = all_trials{1, trial_num}; % channel x samples

    filtered_trial = zeros(num_channels, num_samples);

    for ch = 1:num_channels

        filtered_trial(ch,:) = filter(b, 1, trial(ch, :)); 

            % each trail is 64 * 512 as given in the data set
                                                  
    end
  filtered_data_final{1, trial_num} = filtered_trial;
end


% ----- class 2 -----

for trial_num = 1:num_trials

    trial = all_trials{2, trial_num}; % channel x samples
    filtered_trial = zeros(num_channels, num_samples);

    for ch = 1:num_channels

        filtered_trial(ch,:) = filter(b, 1, trial(ch,:)); 

            % each trail is 64 * 512 as given in the data set
                                                  
    end
  filtered_data_final{2, trial_num} = filtered_trial;
end

% fprintf('Dimensions of X_bandpass: %d x %d \n', size(filtered_data_final));

end

function [X1,X2] = dataAug(arr,num_trials)

%   Increases the effective size of the training dataset

arguments (Input)
    arr
    num_trials
end

arguments (Output)
    X1
    X2
end

fprintf('\n=== SLIDING WINDOW AUGMENTATION ===\n\n');

% Parameters
window_size = 64; % 0.25 seconds
stride = 64;      % 0.25 seconds overlap
fs = 256;

% Calculate how many windows we get per trial

num_samples = 512; % original trial length
num_windows_per_trial = floor((num_samples - window_size)/stride) + 1; %computes to 8

% Apply sliding window to Class 1

class_1_windows = []; % Store all windows from class 1

class_1_windows_2D = []; % Store as a matrix with 2 dimensions

for trial_num = 1:num_trials
    
    trial = arr{1, trial_num}';
    
    % Extract sliding windows from this trial

    for i = 1:num_windows_per_trial
        start_idx = (i-1)*stride + 1;
        end_idx = start_idx + window_size - 1;
        
        % Get window data (64 x 64)
        window_data = trial(start_idx:end_idx, :);
        
        % Store this window
        class_1_windows = cat(3, class_1_windows, window_data);

        class_1_windows_2D = cat(1, class_1_windows_2D, window_data);
        
    end
   
end

% Apply sliding window to Class 2

class_2_windows = []; % Store all windows from class 2

class_2_windows_2D = []; % Store as a matrix with 2 dimensions

for trial_num = 1:num_trials
    
    trial = arr{2, trial_num}';
    
    % Extract sliding windows from this trial
    for i = 1:num_windows_per_trial
        start_idx = (i-1) * stride + 1;
        end_idx = start_idx + window_size - 1;
        
        % Get window data (64 x 64)
        window_data = trial(start_idx:end_idx, :);
        
        % Store this window
        class_2_windows = cat(3, class_2_windows, window_data);

        class_2_windows_2D = cat(1, class_2_windows_2D, window_data);
        
    end

end

X1 = class_1_windows_2D'; 
X2 = class_2_windows_2D'; 

total_windows_class1 = size(class_1_windows, 3);
fprintf('Original trial: 512 samples (2 seconds)\n');
fprintf('Window size: %d samples (%.2f seconds)\n', window_size, window_size/fs);
fprintf('Windows per trial: %d\n', num_windows_per_trial);
fprintf('Original number of trials: 100\n');
fprintf('Total windows created: %d\n', total_windows_class1);
fprintf('Augmentation factor: %.1fx\n', total_windows_class1/100);
fprintf('Each window shape: %d samples x %d channels\n', size(class_1_windows, 1), size(class_1_windows, 2));
% fprintf('Dimensions of X1_bandpass: %d x %d \n', size(X1));
% fprintf('Dimensions of X2_bandpass: %d x %d \n', size(X2));
end

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
    %A = (inv(W))';              % Pattern matrix A
    X1_CSP = (W(: , 1:2))' * X1_mean_centered;
    X2_CSP = (W(: , 1:2))' * X2_mean_centered;

    % fprintf('Dimensions of X1_CSP: %d x %d \n', size(X1_CSP));
    % fprintf('Dimensions of X2_CSP: %d x %d \n', size(X2_CSP));
    % fprintf('Dimensions of W_CSP: %d x %d \n', size(W));
    % fprintf('Dimensions of λ_CSP: %d x %d \n', size(lambda));
    % fprintf('Dimensions of A: %d x %d \n', size(A));
    
end

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

% fprintf('Dimensions of X1_LDA: %d x %d \n', size(X1_LDA));
% fprintf('Dimensions of X2_LDA: %d x %d \n', size(X2_LDA));
% fprintf('Dimensions of W_LDA: %d x %d \n', size(W));
% fprintf('Dimensions of EigenVector matrix: %d x %d \n', size(V));

end

clc,clearvars

path = 'Add your path here';

arr = FIR(path, 100, 45, 1);

[X1, X2] = dataAug(arr, 100);

[X1_CSP, X2_CSP] = CSP(X1, X2, width(X1));

[X1_LDA, X2_LDA] = LDA(X1_CSP, X2_CSP);

disp(X1_LDA);
disp(X2_LDA);
