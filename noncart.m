%% Non-Cartesian MRI using BART

%% Generate trajectories

% Generate k-space trajectory with 64 radial spokes
traj_rad = bart('traj -r -x512 -y64');

figure, scatter(traj_rad(1, :), traj_rad(2,:), 1), 
title('Radial trajectory'), axis off

%% Generate kspace data
% increase the reconstructed FOV a bit
traj_rad2 = bart('scale 0.6', traj_rad);

% simulate eight-channel k-space data
ksp_sim = bart('phantom -k -s8 -t', traj_rad2);

ksp_rss = bart('rss 8', ksp_sim);

figure, scatter(traj_rad(1, :), traj_rad(2,:), 1, log(ksp_rss(:))), 
colormap(gray), title('Kspace data'), axis off

%% Coil-by-coil inverse gridding reconstruction
% inverse gridding
igrid = bart('nufft -i -t', traj_rad2, ksp_sim);

% channel combination
reco1 = bart('rss 8', igrid);

figure, imshow(abs(reco1), [])
title('Inverse Gridding')

%% ESPIRiT calibration
% reconstruct low-resolution image and transform back to k-space
lowres_img = bart('nufft -i -d30:30:1 -t', traj_rad2, ksp_sim);
lowres_ksp = bart('fft -u 7', lowres_img);

% zeropad to full size
ksp_zerop = bart('resize -c 0 308 1 308', lowres_ksp);

% ESPIRiT calibration
sens = bart('ecalib -t 0.02 -m1', ksp_zerop);

figure, imshow3(squeeze(abs(sens)),[])
title('ESPIRiT Maps')

%% Parallel imaging reconstruction

% non-Cartesian parallel imging
reco2 = bart('pics -S -r0.001 -t', traj_rad2, ksp_sim, sens);


figure, imshow(abs(squeeze([reco1, reco2])), []);
title('Inverse Gridding vs Parallel Imaging')


%% Parallel imaging + compressed sensing reconstruction

reco3 = bart('pics -e -l1 -S -r0.001 -t', traj_rad2, ksp_sim, sens);


figure, imshow(abs(squeeze([reco1, reco2, reco3])), []);
title('Inverse Gridding vs Parallel Imaging vs PICS')

