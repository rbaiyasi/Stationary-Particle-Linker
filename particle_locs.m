function [ locs ] = particle_locs( Data1, varargin )
%particle_locs Performs frame-by-frame SNR enhancement and localization of
%particles by radial symmetry method. Currently, SNR boosting is not
%optional, it happens automatically and cannot be tuned.
% INPUT:    Data1: movie of microscope images
%           varagin: { local_thd , Gauss_width , wide2 }
%               local_thd: true/false for whether to use local threshold.
%                   For a good image local threshold is not necessary
%               Gauss_width: estimated Gaussian standard deviation. This is
%                   used to determine the size of the fitting region.
%               wide2: Not entirely sure what this is. "The wide threshold"
%                   wide and wide2 are different. Together they can decide 
%                   how many neighbors to include. More detail is discussed
%                   in Fig 6 in the paper.
% OUTPUT:   locs: cell array of particle localizations in each frame.
[~,~,T] = size(Data1);
locs = cell(1,T);
for t = 1:T
    im = SNR_booster(Data1(:,:,t));
    locs{t} = particle_identify(im,varargin);
end

