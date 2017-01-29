%Copyright (c) 2015 Washington University 
%Created by: Anish Mitra
%

%Washington University hereby grants to you a non-transferable, non-exclusive, royalty-free, 
%non-commercial, research license to use and copy the computer code that may be downloaded within 
%this site (the â€œSoftwareâ€).  You agree to include this license and the above copyright notice in 
%all copies of the Software.  The Software may not be distributed, shared, or transferred to any third party.  
%This license does not grant any rights or licenses to any other patents, copyrights, or other forms of 
%intellectual property owned or controlled by Washington University.  
%If interested in obtaining a commercial license, please contact Washington University's Office of Technology 
%Management (otm@dom.wustl.edu).

%YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED â€œAS ISâ€, 
%WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES 
%OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, 
%COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE 
%OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES 
%ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, 
%WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH 
%DAMAGES. 




%This function uses parabolic interpolation to find the lag using the
%extremum of the lagged cross correlation.

%The function first uses the
%sign of the cross correlation at zero to decide whether to find a maximum
%or minmum. Next, we look for the global max/min.

%lcc is the empirical lagged covariance curve, lags is a vector with the timepoints in each temporal direction (e.g. -8:2:8 for +/- 8 seconds with a 2 second TR). 
%I have set a boundary condition such that any lag greater than 5 seconds is recorded as a NaN-- this is based on our experience that giant lags tend to be noise. You can relax or abolish this 
%boundary condition if you like.

function [peak_lag,peak_cov]= parabolic_interp(lcc,lags)

peak_lag = [];
peak_cov = [];
index = [];
MAX = 5; %Maximal lag to be considered
zero = find(lags==0); %Index for zero lag


%Local Maximum
if lcc(zero) > 0
    [D I] = max(lcc);
    if abs(lags(I)) > MAX
        peak_lag = NaN;
        peak_cov = NaN;
        return
    end
    index = [I-1 I I+1]; %These are the three x values to be used for parabolic interpolation
end

%Local Minimum
if lcc(zero) < 0
    [D I] = min(lcc);
    if abs(lags(I)) > MAX
        peak_lag = NaN;
        peak_cov = NaN;
        return
    end
    index = [I-1 I I+1]; %These are the three x values to be used for parabolic interpolation
end

if isempty(index) == 1
    peak_lag = NaN;
    peak_cov = NaN;
end

if isempty(peak_lag) == 1
    ypoints = lcc(index);
    xpoints = lags(index);
    X = [(xpoints.^2)' xpoints' ones(3,1)]; %Matrix of x-terms
    constants = X\(ypoints)';
    peak_lag = -.5*(constants(2)/constants(1)); %From the first derivative
    peak_cov = constants(1)*peak_lag^2 + constants(2)*peak_lag + constants(3); %peak_cov = abs(peak_cov);
    peak_lag = TR*peak_lag; %Get the lag in seconds
end

end