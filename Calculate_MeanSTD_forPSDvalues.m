function[MeanValue,StdDevValue] = Calculate_MeanSTD_forPSDvalues(PSDvalues)
% This function is being called from scorematic.m file.
% It calculates the Mean and Std Deviation of sigma, theta and delta PSD
% values from the structure passed into the function. PSDvalues are
% calculated in the function called as 'fft_psd_and_statescore_of_epoch'

% Written by Apurva Turakhia. 
% For Sleep and Memory Labs, University of Michigan
% lines 29-33 added 4-24-06 by Theresa to get mean EMG in the PSDvalues excel output file

% calculate the Sigma Mean and Std Deviation
MeanSigma= PSDvalues.sum.sigma / PSDvalues.nPoints.sigma;
arg1=((PSDvalues.nPoints.sigma * PSDvalues.squaresum.sigma) - PSDvalues.sum.sigma.^2);
arg2=PSDvalues.nPoints.sigma * (PSDvalues.nPoints.sigma-1);
StdDevSigma=sqrt(arg1/arg2);

% Calculate the Theta Mean and Std Deviation
MeanTheta= PSDvalues.sum.theta / PSDvalues.nPoints.theta;
arg1=((PSDvalues.nPoints.theta * PSDvalues.squaresum.theta) - PSDvalues.sum.theta.^2);
arg2=PSDvalues.nPoints.theta * (PSDvalues.nPoints.theta-1);
StdDevTheta=sqrt(arg1/arg2);

% Calculate the Delta Mean and Std Deviation
MeanDelta= PSDvalues.sum.delta / PSDvalues.nPoints.delta;
arg1=((PSDvalues.nPoints.delta * PSDvalues.squaresum.delta) - PSDvalues.sum.delta.^2);
arg2=PSDvalues.nPoints.delta * (PSDvalues.nPoints.delta-1);
StdDevDelta=sqrt(arg1/arg2);

% Calculate the EMG Mean and Std Deviation
MeanEMG= PSDvalues.sum.emg / PSDvalues.nPoints.emg;
arg1=((PSDvalues.nPoints.emg * PSDvalues.squaresum.emg) - PSDvalues.sum.emg.^2);
arg2=PSDvalues.nPoints.emg * (PSDvalues.nPoints.emg-1);
StdDevEMG=sqrt(arg1/arg2);

% Put these values in a STRUCT for using it in main program
MeanValue = struct( 'sigma',MeanSigma, 'theta',MeanTheta, 'emg',MeanEMG, 'delta',MeanDelta);
StdDevValue = struct( 'sigma',StdDevSigma, 'theta',StdDevTheta, 'emg',StdDevEMG, ...
    'delta',StdDevDelta);