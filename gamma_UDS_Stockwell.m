%% gamma_UDS_Stockwell
% computes pared-down (but easily extendable) Stockwell transform and
% extracts average power in gamma frequency band

%Requires 'st.m' by Robert Stockwell
%If used, please cite: "Localization of the Complex Spectrum: The S Transform"
% from IEEE Transactions on Signal Processing, vol. 44., number 4, April 1996, pages 998-1001.
%

function [AllMaxGamma] = gamma_UDS_Stockwell(data,SamplingRate)

M=data;
GammaRange = [30 50];
% easy to extend to other frequency bands, e.g.
%     ThetaRange = [3 10];
%     DeltaRange = [0.5 1];


SamplingRateReduction = 1; %no reduction; data are already downsampled
WindowSize = 5;
WindowStep = WindowSize / 2;
TemporalSmoothingTime = 0.5; %
SmoothingKernelPrecision = 1; %sigmas
MaxSamplesPerTrace = size(M, 2);
nEpochsPerTrace=size(M,1);


SamplingRateReduced = SamplingRate / SamplingRateReduction;
params.Fs = SamplingRateReduced;
params.fpass = [0 200];


MaxTraceDuration = MaxSamplesPerTrace/SamplingRate;

nChans=1;

ThisTrace=M; %each epoch is only one channel anyway
ThisTrace = ThisTrace(~isnan(ThisTrace));   %only take real data
nSamplesThisTrace = length(ThisTrace);
ThisTraceDuration = nSamplesThisTrace/SamplingRate;
SamplesThisTrace = length(ThisTrace);
desiredmaxfreq=250;
StockwellMinFreq = 1/(SamplingRateReduced/SamplesThisTrace );
%StockwellMinFreq=1;
StockwellMaxFreq = desiredmaxfreq/(SamplingRateReduced/SamplesThisTrace );
%     StockwellMaxFreq=100;
StockwellFreqStep = 1;

[StockwellSpectro, StockwellTimes, StockwellFreqs] = st(ThisTrace, StockwellMinFreq, ...
    StockwellMaxFreq, 1/SamplingRateReduced, round(StockwellFreqStep / (SamplingRateReduced/SamplesThisTrace)));

% smooth Stockwell spectrum in time domain to reduce salience of short high frequency bouts
TimeStep = median(diff(StockwellTimes));
SmoothingKernelWidth = SmoothingKernelPrecision * TemporalSmoothingTime;
%     TemporalSmoothing = pdf('norm',-SmoothingKernelWidth : TimeStep : SmoothingKernelWidth, 0, TemporalSmoothingTime);
%     TemporalSmoothing = TemporalSmoothing / sum(TemporalSmoothing);

StockwellSpectro = single(abs(st(ThisTrace,StockwellMinFreq,StockwellMaxFreq,...
    1/SamplingRateReduced, round(StockwellFreqStep / (SamplingRateReduced/SamplesThisTrace)) )));
% find max freq in stockwell spectrogram

GammaFreqs = StockwellFreqs >= GammaRange(1) & StockwellFreqs <= GammaRange(2);
IsMaxPowerGamma = bsxfun(@eq, squeeze(StockwellSpectro(GammaFreqs, : )), max(squeeze(StockwellSpectro(GammaFreqs, : )), [], 1));

MaxFreqGamma = sum(bsxfun(@times, StockwellFreqs(GammaFreqs)', IsMaxPowerGamma), 1) ./ sum(IsMaxPowerGamma, 1);

% detect max freq at lower edge of range

LowEdgeGamma = IsMaxPowerGamma(1, :) == true;


AllStockwellSpectro = NaN(length(StockwellFreqs), length(StockwellTimes)*nChans );
AllStockwellTimes = NaN(length(StockwellTimes) *1 , 1);

%on all sweeps, update arrays appropriately
if(size(StockwellSpectro,1)) ~= size(AllStockwellSpectro,1)
    StockwellFreqs
    whos AllStockwellSpectro StockwellSpectro
end

AllStockwellSpectro = StockwellSpectro; %CH, simplified to ONE matrix for testing purposes
AllStockwellTimes = StockwellTimes;
%remove all NaN (likely to be present at end of arrays)
StSpecCleanVect = AllStockwellSpectro(~isnan(AllStockwellSpectro));
nfreqs = size(AllStockwellSpectro,1);
AllStockwellSpectro = reshape(StSpecCleanVect, nfreqs, length(StSpecCleanVect)/nfreqs);
AllStockwellTimes = AllStockwellTimes(~isnan(AllStockwellTimes));

%possible extension: compute max gamma frequencies
GammaFreqs = StockwellFreqs >= GammaRange(1) & floor(StockwellFreqs) <= GammaRange(2);
IsMaxPowerGamma = bsxfun(@eq, squeeze(AllStockwellSpectro(GammaFreqs, : )), max(squeeze(AllStockwellSpectro(GammaFreqs, : )), [], 1));
MaxFreqGamma = sum(bsxfun(@times, StockwellFreqs(GammaFreqs)', IsMaxPowerGamma), 1) ./ sum(IsMaxPowerGamma, 1);
MeanGammaPower = squeeze(mean(AllStockwellSpectro(GammaFreqs, : ), 1));



%optional: compute peak index to quantify 'peakiness' of theta (still needs
%to be optimized)
IsMaxPowerGamma = bsxfun(@eq, squeeze(AllStockwellSpectro(GammaFreqs, : )), max(squeeze(AllStockwellSpectro(GammaFreqs, : )), [], 1));
MaxFreqGamma = sum(bsxfun(@times, StockwellFreqs(GammaFreqs)', IsMaxPowerGamma), 1) ./ sum(IsMaxPowerGamma, 1);

StockwellSpectro=StockwellSpectro';

AllMaxGamma=AllStockwellSpectro([GammaRange(1):GammaRange(2)],:);
AllMaxGamma=AllMaxGamma(IsMaxPowerGamma);
end