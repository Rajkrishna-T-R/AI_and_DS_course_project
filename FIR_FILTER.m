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

fprintf('Dimensions of X_bandpass: %d x %d \n', size(filtered_data_final));

end
