function out = cap_modulator(in, UpSamplingFactor, inphase_pulse,quadrature_pulse)

% real and imaginary separation
inphase_data    = real(in);
quadrature_data = imag(in);

% Upsampling the signal
inphase_data_up     = upsample(inphase_data,UpSamplingFactor);
quadrature_data_up  = upsample(quadrature_data,UpSamplingFactor);

% Filter the data with TX Filter
% inphase_data_up(end+1:end+FilterTaps/2)     = 0;    %zero padding
% quadrature_data_up(end+1:end+FilterTaps/2)  = 0;    %zero padding

inphaseFil      = filter(inphase_pulse,1,inphase_data_up);
quadratureFil   = filter(quadrature_pulse,1,quadrature_data_up);

% inphaseFil      = inphaseFil(FilterTaps/2+1:end);       % eliminate delay
% quadratureFil   = quadratureFil(FilterTaps/2+1:end);    % eliminate delay

% Summation
out = inphaseFil-quadratureFil;


