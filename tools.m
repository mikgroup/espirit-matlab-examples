%% Example 4: Basic Tools
%
% Various tools are demonstrated by manipulating an image.
%
%%

knee_l1 = readcfl('data/knee_l1');

% Zero pad

knee2 = bart('resize -c 1 300 2 300', knee_l1);

% Switch dimensions 1 and 2

tmp = bart('transpose 1 2', knee2);

% Scale by a factor of 0.5

tmp2 = bart('scale 0.5', tmp);

% Join original and the transposed and scaled version along dimension 2.

joined = bart('join 2', knee2, tmp2);

% Flip 1st and 2nd dimension (2^1 + 2^2 = 6)

tmp = bart('flip 6', joined);

% Join flipped and original along dimension 1.

big = bart('join 1', joined, tmp);

% Extract sub-array

tmp = bart('extract 1 150 449', big);
small = bart('extract 2 150 449', tmp);

% Circular shift by 115 pixels

tmp = bart('circshift 1 150', small);
shift = bart('circshift 2 150', tmp);

% Show the final result.

figure, imshow(abs(squeeze(shift)), []);


