clear;
close all;
clc;

%% setting parameters and initializing matrices
tic
TR         = 2;             % TR value 
NP         = 236;           % number of time points
N          = 1000;          % number  of simulated random time series
fShiftStep = 0.1;           % amount of fractional shift
shft       = 0:fShiftStep:5;    % vector containing shift values
MaxLag     = 8;             % Upper bound for lag value in computing lagged correlation function.
dt         = 1;

ESlagMat = zeros(2, size(shft, 2), N);  % for each time series generated, estimated and real shifts are stored in this matrix
SimTsMat = randn(NP, N);                % each column of this matrix contains one simulated random time series

%%%% columns of this matrix will contain the estimated shift for each time series,
%%%% which is then used to calculate mean estimated shift over all randomly generated time series 
AvgEsLag = zeros(N, size(shft, 2));     

%% Introducing a fractional shift and estimating the shift for randomly generated time series

%%%% for each generated time series, a range of fractional shifts ("shft" vector contains 
%%%% the shift values) are introduced using
%%%% fourier transform and the estimated shift is stored in "ESlagMat"
for n = 1:N
    ts      = SimTsMat(:, n);
    
    for st = 1:size(shft, 2)
        
        tmpshft = shft(st);
        
        %%%% Using fourier transform "time shift" property, fractioal shift
        %%%% is introduced in frequency domain. Then using inverse fourier
        %%%% transform, the shifted signal is reconstructed.
        fAxis   = ifftshift((0:NP-1) -ceil((NP-1)/2))/NP/dt;
        ts2     = ifft(fft(ts').*exp(-1i*2*pi*tmpshft*fAxis));
        ts2     = real(ts2);
        s       = [ts, ts2'];
    
        [TD, PC, lagProjection, TDz] = td_mat(s,MaxLag, TR);

        %%%% ESlagMat is a 2*(number of shifts)*(number of simulated time
        %%%% series) matrix. Estimated shift for each time series is stored in a 2*(number of shift values) matrix.
        %%%% the first element will be the real shift and the second
        %%%% element will be the estimated shift
        ESlagMat(1, st, n) = shft(st);
        ESlagMat(2, st, n) = TD(1, 2);
        
        %%%% estimated shifts for time series i are stored in column i of this matrix.
        %%%% then calculating mean value of the estimated shifts is done using a simple mean command
        %%%% which produce a vector containing columnar mean values
        AvgEsLag(n, st)    = TD(1, 2);
        
    end
    
end

%%%% calculating the mean estimated shifts corresponding to a particular
%%%% introduced shift over all time series
AVG = mean(AvgEsLag);

%% visualizing mean estimated shift vs. actual shift

%%%% plot the estimated lag vs the real shift value
%%%% data points are average estimated shifts over radomly generated time
%%%% series. minimum and maximum values of the estimated shifts are
%%%% considered to plot the boundary

figure(1);
plot(shft, AVG/TR, 'c*');
hold on
MIN = min(AvgEsLag); 
MAX = max(AvgEsLag);
plot(shft, MIN/TR, 'b--');
hold on
plot(shft, MAX/TR, 'b--');
xlim([-0.5 5.5]);
ylim([-0.5 5.5]);
axis square
grid on
grid minor
title('estimated vs. actual shift')
xlabel('real shift value')
ylabel('estimated shift value')
hold on

%%%% in the ideal case, points must lay on the y = x line.
%%%% here y = x is plotted in order to have a visual understanding about how
%%%% estimated shift differs from the actual one.
hold on
x = 0:fShiftStep:5;
y = 0:fShiftStep:5;

%%%% here x values are multiplied by TR
plot(y, x, 'k')


toc