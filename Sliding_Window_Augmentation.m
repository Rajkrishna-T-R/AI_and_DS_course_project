function [X1,X2] = dataAug(path,num_trials)
%   Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    path
    num_trials
end

arguments (Output)
    X1
    X2
end

fprintf('\n\n=== SLIDING WINDOW AUGMENTATION ===\n\n');

% Load your dataset
data = load(path);
all_trials = data.eeg_data_wrt_task_rep_no_eog_256Hz_end_trial;

% Parameters
window_size = 64; % 0.25 seconds
stride = 64;       % 0.25 seconds overlap
fs = 256;

% Calculate how many windows we get per trial
num_samples = 512; % original trial length
num_windows_per_trial = floor((num_samples - window_size)/stride) + 1; %computes to 8

%% Apply sliding window to Class 1

class_1_windows = []; % Store all windows from class 1

class_1_windows_2D = []; % Store as a matrix with 2 dimensions

for trial_num = 1:num_trials
    % Get one trial from class 1
    trial = all_trials{1, trial_num}'; % Transpose to 512x64 (time x channels)
    
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

%% Apply sliding window to Class 2

class_2_windows = []; % Store all windows from class 2

class_2_windows_2D = []; % Store as a matrix with 2 dimensions

for trial_num = 1:num_trials
    % Get one trial from class 2
    trial = all_trials{2, trial_num}'; % Transpose to 512x64 (time x channels)
    
    % Extract sliding windows from this trial
    for i = 1:num_windows_per_trial
        start_idx = (i-1)*stride + 1;
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

%% Results
total_windows_class1 = size(class_1_windows, 3);
fprintf('Original trial: 512 samples (2 seconds)\n');
fprintf('Window size: %d samples (%.2f seconds)\n', window_size, window_size/fs);
fprintf('Windows per trial: %d\n', num_windows_per_trial);
fprintf('Original number of trials: 100\n');
fprintf('Total windows created: %d\n', total_windows_class1);
fprintf('Augmentation factor: %.1fx\n', total_windows_class1/100);
fprintf('Each window shape: %d samples x %d channels\n', size(class_1_windows, 1), size(class_1_windows, 2));
fprintf('Dimensions of X1_bandpass: %d x %d \n', size(X1));
fprintf('Dimensions of X2_bandpass: %d x %d \n', size(X2));
end
