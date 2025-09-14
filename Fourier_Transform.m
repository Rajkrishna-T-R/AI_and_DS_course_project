function [] = Fourier_Transform(path,class,trial,channel_req)
% This function dones fourier transform of the data. 
% The function returns the half-spectrum for specified channel.

% Make sure you dont have any file named "plot.m" your current directory
% This function (Fourier Transform) will try to open that if there is file 
% like that which is not what we wanted it to do.
% The "Plot" function in the directory can overide inbuilt plot function


arguments (Input)
    path,
    class,
    trial,
    channel_req
end



data=load(path);

all_trials = data.eeg_data_wrt_task_rep_no_eog_256Hz_end_trial;




num_samples = 512; % ie time duration =2 Sec

Sampling_frq=256; 

Req_signal=all_trials{class,trial}(channel_req,:);

% Example: trail1 = [64 x 512] matrix, where
% - rows: channels
% - columns: samples (time points)



Y=fft(Req_signal,[],2);
% FFT along 2nd dimension (512 samples), []--> automatic size detection


f = (0:num_samples-1)*(Sampling_frq/num_samples);
% Created the frequency vector



  Half_num_samples=floor(num_samples/2)+1;
  Y_half = Y(:, 1:Half_num_samples);
  f_half = f(1:Half_num_samples);
  

  
 
% --- Plotting ---
figure('Color','black');
title(sprintf('Channel %d, Class %d, Trial %d',channel_req,class,trial));
plot(f_half, abs(Y_half), 'y', 'LineWidth', 1.5);
xlabel("Frequency (Hz)", 'Color', 'w');
ylabel("Magnitude", 'Color', 'w');

grid on;

   
    

end
