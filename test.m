% noise_params.mu=50;
% noise_params.mu0=0;
% noise_params.sigma=0;
% noise_params.sigma0=0;
% tic
% [Fsim,~]  = scan_volume(vol_out, PSF_struct, neur_act, ...
%                        scan_params, noise_params, spike_opts, tpm_params); % Perform the scanning simulation
% write_TPM_movie(Fsim, [saveDir,['50_0_0_0.tif']]);
% fprintf('Simulated scanning in %f seconds.\n', toc); 


clean_img=imread('test_clean_160.png');
clean_img=im2double(clean_img);
clean=clean(:);
noise_img=imread('50_0_0_0_160.png');
noise_img=im2double(noise_img);
noise=noise_img-clean_img;
noise=noise(:);
snr(clean,noise)

