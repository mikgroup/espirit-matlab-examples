%%
% non-Cartesian MRI using BART
%%


% Generate k-space trajectory with 64 radial spokes
traj_rad = bart('traj -r -x512 -y64');

% 2x oversampling
traj_rad2 = bart('scale 0.5', traj_rad);

% simulate eight-channel k-space data
ksp_sim = bart('phantom -k -s8 -t', traj_rad2);

% increase the reconstructed FOV a bit
traj_rad2 = bart('scale 0.6', traj_rad);


% inverse gridding
igrid = bart('nufft -i -t', traj_rad2, ksp_sim);

% channel combination
reco1 = bart('rss 8', igrid);


% reconstruct low-resolution image and transform back to k-space
lowres_img = bart('nufft -i -d24:24:1 -t', traj_rad2, ksp_sim);
lowres_ksp = bart('fft -u 7', lowres_img);

% zeropad to full size
ksp_zerop = bart('resize -c 0 308 1 308', lowres_ksp);

% ESPIRiT calibration
sens = bart('ecalib -m1', ksp_zerop);

% non-Cartesian parallel imging
reco2 = bart('pics -S -r0.001 -t', traj_rad2, ksp_sim, sens);



figure, imshow(abs(squeeze([reco1, reco2])), []);
title('Inverse Gridding vs Parallel Imaging')


