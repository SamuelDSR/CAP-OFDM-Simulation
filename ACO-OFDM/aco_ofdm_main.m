clc
clear all;
close all;
addpath('../utils');
disp('-------------ACO-OFDM Simulation-------------------------');
disp('-------------Setting simulation parameters---------------');
tic

%% Parameters
% Using FEC or not
FECcoding = 0;

% Power and Noise
N0          = 10^-21;            % spectral density level of gaussian noise at the receiver
B           = 200*10^6;          % signal bandwidth B = 20MHz
noiseVar    = B*N0;
%AvgPower    = logspace(-3,-1,30); % total available electrical power = E[x(n)_elec^2] or = E[x(k).*conj(x(k))]
AvgPower    = logspace(-3,0,45);
%AvgPower    = 0.01*4;

% DCO-OFDM params
fftSize     = 512;               % fft size
cpSize      = 16;                % cyle prefix size
nSubcar     = fftSize/4;         % number of subcarrier

% LED params
led.minCurrent = 0.1;            % minimun turn-on forward current of LED
led.maxCurrent = 1;              % maximum forward current of LED
led.dcBias     = 0.1;            % DC bias
led.polynomial = [-0.0003,2.0565,-1.0886,0.2855];   % polynomial shaping, (wiener non-linear model)


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
vlcFilterCoeff = 0.4*10^-4; %Attenuated Dirac channel

%% PD filter
%Responsitivity of PD (PD filter considered as a Dirac channel)
pd = 1;

%% Total equivalent channel
totalChannelCoeff        = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%test
%totalChannelCoeff       = ledFilterCoeff;
%totalChannelCoeff       = vlcFilterCoeff;

%Frequency response of total equivalent channel
totalChannelFR           = fft(totalChannelCoeff,fftSize);
if size(totalChannelFR,1) == 1
    totalChannelFR = totalChannelFR';
end

channelStateInformation = totalChannelFR(2:2:2*nSubcar); %Not sure, correct

%Plot Total Channel Frequency Response
% figure
% fvals = B*(-fftSize/2:fftSize/2-1)/fftSize;
% plot(fvals/10^6,fftshift(abs(totalChannelFR)));
% title('Total Channel Frequency response');
% xlabel('Frequency (MHz)');
% ylabel('Magnitude');
% grid on;

%% Bit and power loading (Bit and power loading algorithmes doesn't include the clipping effect)
%================================================================================================================================%
% No bit and power loading
nBitsPerQAM             = 4;
bitAlloc                = nBitsPerQAM*ones(nSubcar,1);
%powerAlloc              = AvgPower*ones(nSubcar,1);
dataRate = sum(bitAlloc)/nSubcar*B*nSubcar/(fftSize+cpSize)
%================================================================================================================================%
% Waterfilling algorithms, only power loading included
% [Capacity, powerAlloc]  = ofdm_waterfilling(nSubcar,totalPower,channelStateInformation,B,N0);
% nBitsPerQAM             = 4;                                              % number of bit per QAM symbol
% %bitAlloc               = nBitsPerQAM*ones(nSubcar,1);
% bitAlloc(powerAlloc==0) = 0;                                              % don't load bits on null power subcarrier

%================================================================================================================================%
% Campello-levin algorithmes
% maxOrder        = 12;
% bitAlloc        = randi([0 maxOrder],nSubcar,1);
% targetBER       = 0.01;
% %gap             = 2/3*erfcinv(targetBER)^2;
% noiseVar        = B*N0;
% gainToNoise     = channelStateInformation.*conj(channelStateInformation)/noiseVar;
% [bitAlloc, powerAlloc] = campello_algo(bitAlloc,targetBER,gainToNoise,nSubcar,4*  AvgPower*nSubcar,maxOrder);


%% Simulation results pre-allocation
snr_elec_tx      = zeros(size(AvgPower));        % Electrical SNR including DC-bias (At the Tx)
snr_elec_nodc_tx = zeros(size(AvgPower));        % Electrical SNR excluding DC-bias (At the Tx)
snr_elec_opt_tx  = zeros(size(AvgPower));        % Optical SNR including DC-bias    (At the Tx)
%snr_elec_rx      = zeros(size(AvgPower));        % Electrical SNR at Receiver

RMSEVM = zeros(nSubcar,length(AvgPower));
numErrors   = zeros(size(AvgPower));
ber         = zeros(size(AvgPower));        % BER
PAPR        = zeros(size(AvgPower));        % PAPR, in dB

%% Simulation Routine
disp('-------------Simulation Routine Started------------------')
nBlockperSimu           = 20000;                                      % number of blocks per simulation unit
nBitsperBlock           = sum(bitAlloc);                            % number of bits per block
nBitsperSimu            = nBitsperBlock*nBlockperSimu;
ErrorTarget             = 10000;                                    % target number of error
maxiBits                = 2*10^8;
disp(strcat('-------------ErrorTarget is:', num2str(ErrorTarget)))
disp(strcat('-------------Maximum simu bits is:', num2str(maxiBits)))

for k = 1:length(AvgPower)
    disp(strcat('-------------AvgPower:', num2str(AvgPower(k)), '------------------'))
    powerAlloc      = AvgPower(k)*ones(nSubcar,1);
    simulations     = 0; % simulation times
    
    while numErrors(k)<=ErrorTarget(k) && simulations*nBitsperSimu <= maxiBits 
        disp(strcat('-------------', num2str(simulations+1), ' Simulation------------------'))
        
        % Binary Input generator
        inData =  random_bin_generator(nBlockperSimu*sum(bitAlloc));

        % Encoding (interleaving / Forward Error Encoding)
        %FEC coding
        if FECcoding ~= 0
            %disp('-------------Forward Encoding----------------------------')
            [encData, CodeTrellis] = conv_encoder(inData,[]);
            
            % interleaving
            tableSize = 10;
            interleaveTable = interleav_matrix(ones(1, tableSize));
            [encData, gain] = interleaving(encData, interleaveTable);
        else
            encData = inData;
        end

        % bit and power loading QAM modulator
        %disp('-------------QAM Modulation and Bit-Power Loading--------')
        modulator = @qam_modulator;
        [qamodData, remBits] = bit_power_loading(encData, bitAlloc, powerAlloc, modulator);
        %mod_signal_var = sum(abs(modData).^2)/length(modData)

        % ACO-OFDM modulator
        %disp('-------------ACO-OFDM Modulation-------------------------')
        acoModParams = [nSubcar,cpSize];
        [modData, blkSize] = aco_ofdm_modulator(qamodData,acoModParams);

         % Electrical Signal Measuring

        elec_signal_mean = sum(modData)/length(modData);
        elec_signal_var  = sum(modData.^2)/length(modData);
        elec_signal_peak = max(abs(modData));
        PAPR(k) = (PAPR(k)*simulations + elec_signal_peak/sqrt(elec_signal_var))/(simulations+1);
        
        % LED channel filtering
        txData = filter(ledFilterCoeff,1,modData);
        
        %LED Clipping and non-linearity distortion
        %disp('-------------LED Clipping------------------')
        %lamda = 1.5;
        %led.dc_bias = lamda*sqrt(P_elec)+led.min %set the dc bias according to the signal variance
        txData  = led_shaping(txData,led);
        %txData = modData;
        
        % Optical Signal Measuring  
        opt_signal_var      = sum(txData.^2)/length(txData);                %optical signal elec after clipping including dc-bias
        opt_signal_nodc     = sum((txData-led.dcBias).^2)/length(txData);   %optical signal elec after clipping excluding dc-bias
        opt_signal_mean     = sum(txData)/length(txData);
        snr_elec_tx(k)      = (snr_elec_tx(k)*simulations + opt_signal_var/noiseVar)/(simulations+1);              
        snr_elec_nodc_tx(k) = (snr_elec_nodc_tx(k)*simulations + opt_signal_nodc/noiseVar)/(simulations+1);         
        snr_elec_opt_tx(k)  = (snr_elec_opt_tx(k)*simulations + opt_signal_mean/noiseVar)/(simulations+1); 
        
        % Led and vlc channel filtering
        %disp('-------------Passing through Channel-------------------')
        rxData              = filter(totalChannelCoeff,1,txData);
        %rxData              = txData;
        
        % Add gaussian noise
        %disp('-------------Adding AWGN Noise-------------------------')
        rxData              = rxData + sqrt(noiseVar)*randn(size(rxData));
        
        % DCO-Demodulation
        %disp('-------------ACO-OFDM Demodulation---------------------')
        demodData           = aco_ofdm_demodulator(rxData,acoModParams);
  
        %Frequency Domaine Equalization
        %disp('-------------Equalization------------------------------')
        %one tap zero-forcing equalization
        demodData = reshape(demodData,nSubcar,nBlockperSimu);
        demodData = zero_forcing(demodData,channelStateInformation);
        
        % EVM initial calculation
        qamodData = reshape(qamodData,nSubcar,nBlockperSimu);
        for m = 1:nSubcar
            tmp = qamodData(m,:)-demodData(m,:);
            RMSEVM(m,k) = RMSEVM(m,k)+sum(tmp.*conj(tmp));
        end    
        demodData = reshape(demodData,nSubcar*nBlockperSimu,1);
             
        % QAM Demodulation
        %disp('-------------QAM Demodulation--------------------------')
        demodulator         = @qam_demodulator;
        demodData           = de_bit_power_loading(demodData, bitAlloc, powerAlloc, demodulator);
       
        % Decoding
        if FECcoding ~= 0
            disp('-------------------Decoding----------------------------')
            deleaving = @de_interleaving;
            deleav_data = vlc_decode(demodData, deleaving, interleaveTable);
            
            decoder = @conv_decoder;
            decoderData= vlc_decode(deleav_data, decoder, CodeTrellis);
        else
            decoderData = demodData;
        end

        % Calculation of num of errors
        %disp('-------------BER Calculation-----------------')
        [number_of_errors,~] = biterr(inData,decoderData);
        numErrors(k) = numErrors(k) + number_of_errors;
        simulations = simulations+1;
    end
    
     % EVM calculation fix
    RMSEVM(:,k) = (RMSEVM(:,k)/(nBlockperSimu*simulations)./powerAlloc).^0.5;
    
    snr_elec_tx(k)           = 10*log10(snr_elec_tx(k));
    %snr_elec_rx(k)           = 10*log10(snr_elec_rx(k));
    snr_elec_nodc_tx(k)      = 10*log10(snr_elec_nodc_tx(k));
    snr_elec_opt_tx(k)       = 10*log10(snr_elec_opt_tx(k));

    PAPR(k)             = 10*log10(PAPR(k));
    ber(k)              = numErrors(k)/(simulations*nBitsperBlock*nBlockperSimu);
    
end
toc

%% thereoctical BER
% EbNo = snrr - 10*log10(bitAlloc(1));
% ber_t = berawgn(EbNo,'qam',2^bitAlloc(1))

%% Plot
figure;
semilogy(snr_elec_tx,ber,'r');

figure;
semilogy(snr_elec_nodc_tx,ber,'g');

figure;
semilogy(snr_elec_opt_tx,ber,'b');

figure
semilogy(AvgPower, theory_ber_rectqam(2^nBitsPerQAM,1./RMSEVM(1,:).^2)*nBitsPerQAM,'r');
hold on
semilogy(AvgPower, ber,'b');