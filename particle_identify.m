function params = particle_identify(im, varargin)
% particle_identify Identifies particles using SNR boosted image data.
% Original code by Dr. Bo Shuang. Updated by Rashad Baiyasi
% INPUT:    im: single frame array, should have been SRN boosted previously
%           varagin: { local_thd , Gauss_width , wide2 }
%               local_thd: true/false for whether to use local threshold.
%                   For a good image local threshold is not necessary
%               Gauss_width: estimated Gaussian standard deviation. This is
%                   used to determine the size of the fitting region.
%               wide2: Not entirely sure what this is. "The wide threshold
%                   wide and wide2 are different. Together they can decide 
%                   how many neighbors to include. More detail is discussed
%                   in Fig 6 in the paper."
%% initial conditions & input parameters
defargs = { true , 3 , 2 }; % { local_thd , Gauss_width , wide2 }
if ~isempty(varargin)
    arginds = find(~cellfun(@isempty,varargin));
    defargs(arginds) = varargin(arginds);
end
[local_thd, Gauss_width, wide2] = defargs{:};
wide = floor(wide2);
n = 3; % how many std to add up as a threshold
% FWHM = 2.35*Gauss_width

%% calculate threshold map
w = 50;% the local region to calculate the local background and threshold
% usually 50X50 and shift by 25 is a good choice for a 512X512 image
[Y,X] = size(im);
counts = zeros(Y, X);% counts store the times each pixel has contributed in
% local background calculation
bg = counts; % record the local background
sd = counts; % record the local standard deviation
if local_thd == true
for i = 1 : ceil(Y / w * 2) - 1
    for j = 1 : ceil(X / w * 2) - 1
        % 1, select the local region
        % 2, sort the pixels based on their intensities
        % 3, find the 50% and 75% point as discussed in the paper
        % 4, calculate local background(bg), standard deviation(sd) and
        % count
        
    % edited by Rashad Baiyasi 4/3/17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpRows = 1+w/2*(i-1):min(Y,w/2*(i+1)); %added to reduce computation 
        tmpCols = 1+w/2*(j-1):min(X,w/2*(j+1)); %time and improve readability
        im_local = im(tmpRows, tmpCols);
        im_local = sort(im_local(:));
        n_loc = numel(im_local);
        bg(tmpRows, tmpCols) = bg(tmpRows, tmpCols) ...
            + im_local(round(n_loc/2));
        %sd = 0.82 to 0.5 of cumulative distribution
        sd(tmpRows, tmpCols) = sd(tmpRows, tmpCols) ...
            + im_local(round(n_loc*0.5)) - im_local(round(n_loc*0.18));
        counts(tmpRows, tmpCols) = counts(tmpRows, tmpCols) + 1;
    % end edits 4/3/17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end % for j
end % for i
bg = bg ./ counts;
sd = sd ./ counts;
thd_map = bg + n * sd;% determine the local threshold
else
    % calculate the background and standard deviation once for the whole
    % image
    im_local = sort(im(:));
    n_loc = numel(im_local);
    bg(1:Y, 1:X) = im_local(round(n_loc/2));
    sd(1:Y, 1:X) = im_local(round(n_loc/2)) - im_local(round(n_loc/4));
    thd_map = bg + n * sd;
end % if local_thd

%% calculate max_map
% the idea of this part also originates from Arnauld Serge's work.
center = im(1+wide:Y-wide, 1+wide:X-wide); %remove image borders
max_map = zeros(Y, X); %size of full image
pos_check = max_map;
if numel(thd_map) > 1
    thd_map = thd_map(1+wide:Y-wide, 1+wide:X-wide);
end
max_map(1+wide:Y-wide, 1+wide:X-wide) = center > thd_map;% each center is compared
% with its local threshold.
% the two loops below are intended to select the local maximums that meet two
% conditions: 
% 1, the selected neighbors are not brighter than the center;
% 2, the selected neighbors are also brighter than the local threshold.
%RB - ???? what is going on here?
for i = -wide : wide
    for j = -wide : wide
        if i^2 + j^2 <= wide2^2
            %RB - not sure what this is doing... why +i to beginning & end?
            pos_check(1+wide:Y-wide, 1+wide:X-wide) = ...
                im(1+wide+i:Y-wide+i, 1+wide+j:X-wide+j) <= center & ...
                im(1+wide+i:Y-wide+i, 1+wide+j:X-wide+j) > thd_map;
            max_map = max_map .* pos_check;
            pos_check = zeros(Y, X);
        end
    end
end
max_map = max_map .* (im - bg);
%
%% calculate the subpixel position of each particle
% We use Parthasarathy's radial symmetry method here because it is fast and
% accurate. Detail see nmeth.2071
match_r = 2 * Gauss_width; % the size of fitting region is match_r*2+1
% disp(['match_r = ',num2str(match_r)]);
if isnan(sum(max_map(:))) || isinf(sum(max_map(:)))
    params = [];
    return
end
params = zeros(ceil(sum(max_map(:))),3);
k = 1;
sig_thd = 0.204*(2*match_r+1)*0.9;% threshold based on the width of noise. 
% Width of a real particle must smaller than 90% of the width of noise.
% the magic number 0.204 is the average value per pixel averaged by 1000
% trails
while sum(max_map(:)~=0) > 0
    [~, q] = max(max_map(:));
    q = q(1);
    max_map(q) = 0;
    row = ceil(q / Y);
    row1 = max(row - match_r, 1);
    row2 = min(row + match_r, X);
    clm = q - (row - 1) * Y;
    clm1 = max(clm - match_r, 1);
    clm2 = min(clm + match_r, Y);
    [xc yc sigma] = radialcenter(im(clm1:clm2,row1:row2));
    if sigma < sig_thd
        params(k, 1:3) = [xc+row1-1, yc+clm1-1, sigma];
        k = k + 1;
    end
end
params(k:end,:) = [];
%}