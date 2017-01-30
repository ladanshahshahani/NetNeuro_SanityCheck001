function [ TD, PC , lagProjection, TDz ] = td_mat( Roi_Ts, maxLag, TR)

%   TD_MAT constructs time delay and zero-mean time delay matrices given ROI
%   time series. Time delay matrix is a matrix constructed using parabolic
%   interpolation function(parabolic_interp.m).
%   First it computes the pairwise cross-correlation function for all combinations of ROIs' time series.
%   Then using parabolic interpolation, it computes the peak of cross-correlation function and fills up the TD matrix.
%   Then subtracting the mean from each column of the TD matrix, it constructs TDz matrix. 
%   It also computes and returns the correlation value at the peak.

%   INPUTS:
%   Roi_Ts        : matrix of ROI time series
%   TR            : time of repetition for data acquisition
%   maxLag        : Upper bound for lag value in computing lagged correlation function. 

%   OUTPUTS:
%   TD            : time delay matrix
%   PC            : peak correlation matrix
%   lagProjection : average of each columns of TD matrix
%   TDz           : zero-meaned TD

%% Part1: computing the pairwise cross-correlation function
nnodes    = size(Roi_Ts, 2);
TS        = Roi_Ts;

%%%% applying xcov on time series and calculating cross-correlation
%%%% function (r) with the corresponding lag values (lags)

%%%% [C, Lags] = xcov(A, MAXLAG, 'option') estimates Cross-covariance function for all
%%%% possible combinations of columns of A over the range [-MAXLAG,
%%%% +MAXLAG]. 'biased' option scales the raw cross-covariance by length of the time series. 
[r, lags] = xcov(TS, maxLag, 'biased');

%%%% reshape cross-correlation functions (r) to a 3D matrix to handle the
%%%% function easier
R         = reshape(r, (2*maxLag)+1, nnodes, nnodes);

%%%% initializing time delay matrix (TD) and peak correlation matrix (PC)
TD        = zeros(nnodes);
PC        = zeros(nnodes);

lags = lags';

%%%% zero padd reshaped cross correlation function (R) and lag vector (lags) to avoid indexes equal to zero in the parabolic_interp.m routine
% R    = padarray(R, [1 0 0], 'pre');
% lags = padarray(lags, 1, 'pre');

%% Part2: parabolic interpolation 
for i = 1:nnodes
    for j = 1:nnodes
        [peak_lag,peak_cov]= parabolic_interp_me((R(:, i, j))',lags, TR);               
        
        %%%% column i of TD is a lag map with respect to node i. 
        TD(i, j) = peak_lag;
        PC(i, j) = peak_cov;
    end
end

%%%% mean of TD's columns gives us the lag projection map
lagProjection = nanmean(TD);

%%%% subtracting mean of each column of TD from elements of the corresponding column
%%%% column i of TDz is the zero-centered lag map with respect to node i  
TDz           = bsxfun(@minus,TD,nanmean(TD));

end
