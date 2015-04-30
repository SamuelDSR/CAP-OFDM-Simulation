clear
clc


%% test 1

% s0 = randi([1 16],100,1);
% ls = length(s0);
% h1 = randi([2 9],8,1);
% lh1 = length(h1);
% h2 = randi([2 9],8,1);
% lh2 = length(h2);
% 
% cp = 14;
% cps = zeros(ls+cp,1);
% cps(1:cp) = s0(end-cp+1:end);
% cps(end-ls+1:end) = s0;
% 
% cps(end+1:end+20) =rand(20,1);
% 
% y1 = filter(h1,1,cps);
% 
% y1t = ifft((fft(s0).*fft(h1,ls)),ls);
% 
% % test
%  test0 = y1(cp+1:end-20)-y1t;
% % test1 = ifft((fft(y1(lh1+1:end),ls)./fft(h1,ls)),ls) - s0;
% 
% % y2 = filter(h2,1,y1);
% % 
% % y2t = ifft((fft(y1(cp+1:end),ls).*fft(h2,ls)),ls);
% % 
% % y2(cp+1:end)-y2t
% 
% % test 
% % y2t = filter(conv(h1,h2),1,cps);
% % y2-y2t


%% test 2
FilterSpan          = 8;
UpSamplingFactor    = 4;
FilterTaps          = UpSamplingFactor*FilterSpan;

wt = 0:FilterSpan*2*pi/(FilterTaps):FilterSpan*2*pi;    % time axis for pulse generation (one period)  
rolloff = 0.1;
filtDef = fdesign.pulseshaping(UpSamplingFactor,'Square Root Raised Cosine','N,Beta',FilterTaps,rolloff);
rrcFilter = design(filtDef);
% fvtool(rrcFilter);

% normalize the filter coefficient
inphase_pulse       = rrcFilter.Numerator.* cos(wt);        % in-phase pulse
quadrature_pulse    = rrcFilter.Numerator.* sin(wt);        % quadrature pulse

in = randi([1 9],100,1);
in(end+1:end+FilterSpan) = 0;
in=[in;in];

qd = randi([1 9],100,1);
qd(end+1:end+FilterSpan) = 0;
qd=[qd;qd];

ingain = max(abs(conv(inphase_pulse,inphase_pulse)));
qdgain = max(abs(conv(quadrature_pulse,quadrature_pulse)));

%plot(conv(inphase_pulse,inphase_pulse))

in_s = upsample(in,UpSamplingFactor);
qd_s = upsample(qd,UpSamplingFactor);

tx1 = filter(inphase_pulse,1,in_s);
tx2 = filter(quadrature_pulse,1,qd_s);
tx  = tx1-tx2;

complexfilter = inphase_pulse/ingain-1i*quadrature_pulse/qdgain;
rx1  = ifft(fft(tx).*fft(complexfilter',length(tx)),length(tx));


rx2 = filter(inphase_pulse,1,tx) +1i*filter(quadrature_pulse,1,tx);

res1 = rx1-rx2

indownsample = zeros(100,1);
qddownsample = zeros(100,1);

delay = FilterSpan*UpSamplingFactor;
for i = 1:100
    tmp = rx1(432+delay+1+(i-1)*UpSamplingFactor)
    indownsample(i) = real(tmp);
    qddownsample(i) = imag(tmp);
end

resin=round(indownsample)-in(109:208);
resqd=round(qddownsample)-qd(109:208);
scatter(indownsample,qddownsample)


