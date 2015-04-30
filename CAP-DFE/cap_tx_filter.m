function [inphase_pulse, quadrature_pulse] = cap_tx_filter(FilterSpan,UpSamplingFactor)

FilterTaps          = UpSamplingFactor*FilterSpan;

wt = 0:FilterSpan*2*pi/(FilterTaps):FilterSpan*2*pi;    % time axis for pulse generation (one period)  
rolloff = 0.1;
filtDef = fdesign.pulseshaping(UpSamplingFactor,'Square Root Raised Cosine','N,Beta',FilterTaps,rolloff);
rrcFilter = design(filtDef);
% fvtool(rrcFilter);

% normalize the filter coefficient
inphase_pulse       = rrcFilter.Numerator.* cos(wt);        % in-phase pulse
quadrature_pulse    = rrcFilter.Numerator.* sin(wt);        % quadrature pulse

inphase_pulse       = inphase_pulse/max(abs(inphase_pulse));
quadrature_pulse    = quadrature_pulse/max(abs(quadrature_pulse));

%%plot shaping filter and matched filer
figure
plot(inphase_pulse);
hold on;
plot(quadrature_pulse);

% figure
% plot(inphase_inphase_pulse)
% hold on
% plot(quadrature_quadrature_pulse)
% plot(inphase_quadrature_pulse) 