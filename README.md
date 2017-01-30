# NetNeuroProject
Sanity check lag extraction
"demo_SynthTimeSeries.m" is intended for checking if the artificially introduced lags are correctly estimated.
Both full and fractional delay are introduced. It calls "td_mat.m", which in turn calls "parabolic_interp_me.m".
"td_mat.m" is intended for constructing time delay (TD), peak correlation (PC), zero-centered time delay (TDz), and
mean over columns of TD (LagProjection) matrices. Its inputs are : matrix containing ROI time series (Roi_Ts), repetition time of scanning (TR), 
and the maximum desired lag for computation of lagged cross-correlation function (MaxLag).
In order to estimated the lag between all possible paire of time series, it uses "parabolic_interp_me.m".
"parabolic_interp_me.m" is my modified version of "parabolic_interp.m" code developed by Anish Mitra to calculate lag using
parabolic interpolation. 
To see the results you must run "demo_SynthTimeSeries.m".
"demo_SimSynthTimeSeries.m" is intended for plotting [estimated lag vs. actual lag] to capture the variations of the estimation.
first it generates N random time series. For each of these time series a range of lags is created. You may change upper and lower bounds 
of the range and fractional lag steps.
"parabolic_interp.m" is the original code written by Anish Mitra. I have made some modifications to the code. The modifications are:
1. at line 44 of "parabolic_interp.m" I added zero(1) = [];
2. at line 59 of "parabolic_interp.m", I replaced '[D, I] = min(lcc)' with the following lines:
  lcc2 = -1*lcc;
  lcc3 = lcc2(2:end);
  [D, I] = min(lcc3);
  I = I+1
  the last 3 commands I added are basically for handling the zero index that is added to both "lags" vector and "R" in the "td_mat.m" routine (lines 45, 46)
3. at line 76 of "parabolic_interp.m" transpose signs are omitted.
