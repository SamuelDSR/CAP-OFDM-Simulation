clc
clear all;
close all;
addpath('../utils');
disp('-------------CAP-FDE Simulation-------------------------');
disp('-------------Setting simulation parameters---------------');
tic

%% Set up Parameters
% Power and Noise
N0                  = 10^-21;                   % spectral density level of gaussian noise at the receiver
B                   = 100*10^6;                  % signal bandwidth B = 100MHz 
noise_var           = B*N0;
%AvgPower            = 0.025;                   % total available electrical power = E[x(n)_elec^2] or = E[x(k).*conj(x(k))]
AvgPower    = logspace(-3,-1,30);
%led params
led.minCurrent      = 0.1;                          % minimun turn-on forward current of LED
led.maxCurrent      = 1;                            % maximum forward current of LED
led.dcBias          = .5;                           % DC bias
led.polynomial      = [-0.0003,2.0565,-1.0886,0.2855];   % polynomial shaping, (wiener non-linear model)
%led.polynomial = [0.0771,1.1760];                   % linear model

%% Cap-FDE params
FilterSpan          = 8;
UpSamplingFactor    = 4;
[inphase_pulse, quadrature_pulse] = cap_tx_filter(FilterSpan,UpSamplingFactor);
inphasegain    =  max(abs(conv(inphase_pulse,inphase_pulse)));
quadraturegain =  max(abs(conv(quadrature_pulse,quadrature_pulse)));
complexMatchFilter  = inphase_pulse/inphasegain - 1i*quadrature_pulse/quadraturegain;

%% Block organization
nBitsPerSymbol      = 2;                        %QAM order               
paddingzeros        = FilterSpan;
blockSize           = 256;
nSymbolPerBlock     = blockSize-paddingzeros;   
dataRate = nBitsPerSymbol*B*nSymbolPerBlock/blockSize

%% LED filter (frequency response of LED)
% low pass filter type impulse response h(t) = exp(-2*pi*fc*t) where fc is the 3dB cutoff frequency
% for fc = 5MHz, fs = 100MHz, 32 taps is already enough to simulate this low pass channel (the last filter coefficient is already [[1.58915973878541e-05]])
% 16 taps of CP is enough to composent this channel because filter coefficient at 16 tap is already [[0.00242197535626728]]
ledCutoffFrequecy   = 20*10^6;
ledFilterOrder      = 31;             %
ledFilterCoeff      = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);
% plot(ledFilterCoeff);
% title('Led LPF filter');
% grid on;

%% VLC channel filter
%Dirac channel type
vlcFilterCoeff = 4*10^-5; %Attenuated Dirac channel

%% PD filter
%Responsitivity of PD (PD filter considered as a Dirac channel)
pd = 1;

%% Total equivalent channel
totalChannelCoeff        = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%test
%totalChannelCoeff       = ledFilterCoeff;
%totalChannelCoeff       = vlcFilterCoeff;

%Frequency response of total equivalent channel
fftSize = blockSize*UpSamplingFactor;
totalChannelFR           = fft(totalChannelCoeff,blockSize*UpSamplingFactor);
if size(totalChannelFR,1) == 1
    totalChannelFR = totalChannelFR';
end

%Plot Total Channel Frequency Response
% figure
% fvals = B*(-fftSize/2:fftSize/2-1)/fftSize;
% plot(fvals/10^6,fftshift(abs(totalChannelFR)));
% title('Total Channel Frequency response');
% xlabel('Frequency (MHz)');
% ylabel('Magnitude');
% grid on;

%% Simulation results pre-allocation
snr_elec_tx      = zeros(size(AvgPower));        % Electrical SNR including DC-bias (At the Tx)
snr_elec_nodc_tx = zeros(size(AvgPower));        % Electrical SNR excluding DC-bias (At the Tx)
snr_elec_opt_tx  = zeros(size(AvgPower));        % Optical SNR including DC-bias    (At the Tx)
%snr_elec_rx      = zeros(size(AvgPower));        % Electrical SNR at Receiver

numErrors   = zeros(size(AvgPower));
ber         = zeros(size(AvgPower));        % BER
PAPR        = zeros(size(AvgPower));        % PAPR, in dB

%RMSEVM = zeros(length(AvgPower),1);

%% Simulation Routine
disp('-------------Simulation Routine Started------------------') 
nBlockperSimu           = 20000;                                        % number of blocks per simulation unit 
nBitsperBlock           = nBitsPerSymbol*nSymbolPerBlock;               % number of bits per simulation unit
nBitsperSimu            = nBitsperBlock*nBlockperSimu;
ErrorTarget             = 20000*ones(size(AvgPower));                                   % target number of error
maxiBits                = 2*10^9;
disp(strcat('-------------ErrorTarget is:', num2str(ErrorTarget)))
disp(strcat('-------------Maximum simu bits is:', num2str(maxiBits)))

for k = 1:length(AvgPower)
    disp(strcat('-------------AvgPower:', num2str(AvgPower(k)), '------------------'))
    simulations     = 0; % simulation times

    while numErrors(k)<=ErrorTarget(k) && simulations*nBitsperSimu <= maxiBits
            disp(strcat('-------------', num2str(simulations+1), ' Simulation------------------'))
            % Generate random signal
            inData =  random_bin_generator(nBitsperBlock*nBlockperSimu);
            
            % QAM modulation
            ModOrder        = 2^nBitsPerSymbol;
            qamodData       = qam_modulator(inData, ModOrder);
            
            % Block transmission and adding cycle Prefix
            modData         = cap_fde_block(qamodData,nSymbolPerBlock,paddingzeros);
            
            % CAP modulation
            modData = cap_modulator(modData,UpSamplingFactor,inphase_pulse,quadrature_pulse);
            
            % Electrical Signal Measuring
            elec_signal_mean = sum(modData)/length(modData);
            elec_signal_var  = sum(modData.^2)/length(modData);
            elec_signal_peak = max(abs(modData));
            PAPR(k) = (PAPR(k)*simulations + elec_signal_peak/elec_signal_var)/(simulations+1);
            alpha = sqrt(AvgPower(k)/elec_signal_var);
            modData = alpha*modData;
            
            % Led channel filtering
            txData = filter(ledFilterCoeff,1,modData);
            %rxData  = txData;
            
            % LED Clipping and non-linearity distortion
            txData  = led_shaping(txData, led);
            
            % Optical Signal Measuring
            opt_signal_var      = sum(txData.^2)/length(txData);                %optical signal elec after clipping including dc-bias
            opt_signal_nodc     = sum((txData-led.dcBias).^2)/length(txData);   %optical signal elec after clipping excluding dc-bias
            opt_signal_mean     = sum(txData)/length(txData);
            snr_elec_tx(k)      = (snr_elec_tx(k)*simulations + opt_signal_var/noise_var)/(simulations+1);
            snr_elec_nodc_tx(k) = (snr_elec_nodc_tx(k)*simulations + opt_signal_nodc/noise_var)/(simulations+1);
            snr_elec_opt_tx(k)  = (snr_elec_opt_tx(k)*simulations + opt_signal_mean/noise_var)/(simulations+1);
            
            % vlc channel filtering
            rxData = filter(vlcFilterCoeff,1,txData);
            
            % Add gaussian noise
            detectData = rxData + sqrt(noise_var)*randn(size(rxData));
            %detectData = rxData;
            
            % Deblock transmission
            %deblkData = deblock_transmission(detectData,blockSize*UpSamplingFactor,paddingzeros*UpSamplingFactor);
            
            % Remove dc-bias
            detectData = detectData - vlcFilterCoeff*(1.1760*led.dcBias+0.0771);
            
            % CAP-FDE demodulation
            channelStateInformation = 1.1760*alpha*totalChannelFR;
            equalizedData = cap_fde_demodulator(detectData, nSymbolPerBlock, FilterSpan, UpSamplingFactor, complexMatchFilter, channelStateInformation);
            
            % EVM initial calculation    
            %RMSEVM(k,1) = RMSEVM(k,1)+sum(tmp.*conj(tmp));
            
          %% QAM Demodulation
            demodData = qam_demodulator(equalizedData, ModOrder);
            
            % Calculation of num of errors
            [number_of_errors,~] = biterr(inData,demodData);
            numErrors(k) = numErrors(k) + number_of_errors;
            simulations = simulations+1;
    end
    
    % EVM calculation fix
    %RMSEVM(k,1) = (RMSEVM(k,1)/(nSymbolPerBlock*nBlockperSimu*simulations)).^0.5;
    
    snr_elec_tx(k)           = 10*log10(snr_elec_tx(k));
    %snr_elec_rx(k)           = 10*log10(snr_elec_nodc_tx(k)*channel^2);
    snr_elec_nodc_tx(k)      = 10*log10(snr_elec_nodc_tx(k));
    snr_elec_opt_tx(k)       = 10*log10(snr_elec_opt_tx(k));

    PAPR(k)             = 10*log10(PAPR(k));
    ber(k)              = numErrors(k)/(simulations*nBitsperBlock*nBlockperSimu);
    
end
toc

%% Plot
figure;
semilogy(snr_elec_tx,ber,'r');

figure;
semilogy(snr_elec_nodc_tx,ber,'g');

figure;
semilogy(snr_elec_opt_tx,ber,'b');

figure
%semilogy(AvgPower, theory_ber_rectqam(2^nBitsPerSymbol,1./RMSEVM.^2),'r');
%hold on
semilogy(AvgPower, ber,'b');
