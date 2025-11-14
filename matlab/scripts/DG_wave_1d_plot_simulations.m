clc;clear;close all;

save_video = false;
plot_images = true;

% settings
filename_index = 205;
plot_times = 0:0.1:1;

% read file
filepath = "data/matlab/wave/hdf5/";
filename = filepath + "wave-" + filename_index + ".h5";
wave_postprocessor = dg1d.WavePostprocessor1D();
wave_postprocessor.read_from_hdf5(filename);
disp(wave_postprocessor.metadata)

% plot solution
if plot_images
    wave_postprocessor.plot_solutions(plot_times);
    drawnow;
end

% save video
video_filepath = "data/matlab/wave/videos/";
video_filename = video_filepath + "wave-" + filename_index + ".avi";
if save_video
    wave_postprocessor.create_and_save_wave_animation(video_filename, 2000);
end