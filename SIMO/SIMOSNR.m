% SIMO SNR
% REU Supervisor: Dr. Pan
% REU Mentor: Chenpei Huang
% Written By: Nhat Nguyen

clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates how multiple antennas in a SIMO system could
% improve SNR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTROL PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna
antennaCount = 10; % the number of antennas in this SIMO system

% SNR
channelFadeStd = 1; % channel variance, used to generate channel response
noiseStd = 1; % noise power is this squared
transmitterAveragePower = 20; % transmitter average power (Watts)

% Subcarrier
constellationCount = 16; % the mapping of the constellation, must be a power of 2
bitPerConstellation = log2(constellationCount);

% OFDM
totalSubcarriers = 64; % total number of subcarriers in this OFDM system
pilotSubcarriers = 0; % total number of subcarriers reserved for pilot signals, set to 0 to remove pilot signaling
dataSubcarriers = totalSubcarriers - pilotSubcarriers; % the remaining subcarriers are used for data payload

% Data
symbolsPerFrame = 100; % the number of symbols per frame before cyclic prefix is wrapped around
cyclicLength = 4; % N symbols will be copied and wrapped around the transmitted signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stage 1: Setup
% Generate the Channel Response
% generate a complex gaussian-distributed value for each subcarrier 
% received in each antenna

h = zeros([ totalSubcarriers, antennaCount ]);
for n = 1:totalSubcarriers
    for m = 1:antennaCount
        h(n,m) = h(n,m) + (normrnd(0, channelFadeStd/sqrt(2)) + 1i * normrnd(0, channelFadeStd/sqrt(2))); % cfr 
    end
end


% Generate the Data Transmit Payload
% assume that the binary data is truly random with an even distribution,
% this would give all constellation mappings equally-likely probablies,
% meaning the transmitter's average power is the average of the
% constellation point's magnitude squared
binaryData = randi([0, 1], [ bitPerConstellation * symbolsPerFrame * ( dataSubcarriers ) , 1]);

% Subcarrier Mapping
dataIndices = 1:totalSubcarriers; % create an array of subcarriers
if pilotSubcarriers > 0
    pilotIndices = 1:totalSubcarriers/pilotSubcarriers:totalSubcarriers;
else
    pilotIndices = [];
end

dataIndices(pilotIndices) = []; % remove pilot indices

%% Stage 2: Create Transmit Signal
% QAM Modulation
QAMTXSignal = qammod(binaryData, constellationCount,'InputType','bit','UnitAveragePower',true) * sqrt(transmitterAveragePower); 

% Reshape Data
reshapedTXSignal = reshape(QAMTXSignal, dataSubcarriers, symbolsPerFrame);
mappedTXSignal = zeros([totalSubcarriers, symbolsPerFrame]);

% Mapping Data onto Subcarriers
for i = 1:(dataSubcarriers)
    mappedTXSignal(dataIndices(i),:) = reshapedTXSignal(i,:);
end

% Mapping Pilot onto Subcarriers
for i = 1:(pilotSubcarriers)
    mappedTXSignal(pilotIndices(i),:) = (1+1i) * sqrt(transmitterAveragePower);
end

ifftTXSignal = ifft(mappedTXSignal); % perform IFFT

% Cyclic Prefix (Final Signal)
Txsignal = [ ifftTXSignal((totalSubcarriers - cyclicLength + 1 : end), :); ifftTXSignal ];
%% Stage 3: Received Signal Environment & Decoding
% make four copies of the signal, each for a different antenna, pass it
% through different channel responses, decode

antennaQAM = zeros([ dataSubcarriers, symbolsPerFrame, antennaCount ]);
antennaSNR = zeros([ dataSubcarriers , antennaCount]); % matrix of size M x N where M is the subcarrier index and N is the antenna index

for i = 1:antennaCount
    % pass each signal through its own independent "channel response", this
    % simulates multiple antennas placed at different locations
    % in addition, noise is added
    
    mappedRXSignal = h(:,i) .* mappedTXSignal; % channel response
    for n = 1:totalSubcarriers
        for m = 1:symbolsPerFrame
            mappedRXSignal(n,m) = mappedRXSignal(n,m) + (normrnd(0, noiseStd/sqrt(2)) + 1i * normrnd(0, noiseStd/sqrt(2))); % noise
        end
    end

    %Equalization (Based on Channel estimation)
    mappedRXSignal_eq = mappedRXSignal./h(:,i);
    
    reshapedRXSignal = mappedRXSignal_eq(dataIndices,:); % remove pilot indices
    
    antennaQAM(:,:,i) = reshapedRXSignal;
    QAMRXSignal = reshape(reshapedRXSignal,[],1); % reshape to QAM vector

    %*
    % demod QAM
    binaryRXData = qamdemod(QAMRXSignal / sqrt(transmitterAveragePower), constellationCount,'OutputType','bit','UnitAveragePower',true);
    
    % calculate error and SNR
    antennaSNR(:,i) = 10 * log10(transmitterAveragePower .* abs(h(dataIndices,i)).^2 / noiseStd^2);
    temp = 10 * log10(transmitterAveragePower .* mean(abs(h(dataIndices,i)).^2) / noiseStd^2);
    
    error = sum(abs(binaryData - binaryRXData),'all') / (bitPerConstellation * symbolsPerFrame * ( dataSubcarriers ));
    fprintf("Antenna %d, Bit Error: %.3f%%, Average SNR: %.2f [dB]\n", i, error*100, temp);
end

%% Stage 4: Signal Selection and Decoding Based On Best SNR
% Preallocate
bestQAMSignal = zeros([ dataSubcarriers, symbolsPerFrame ]); % signal is in QAM
equalizer = zeros([ dataSubcarriers, 1]);
[bestCFR, bestAntenna] = max(abs(h(dataIndices,:).^2)');% Best channel frequency response(squared)
bestSNR = 10 * log10(transmitterAveragePower * bestCFR / noiseStd^2);

for num_sc = 1 : dataSubcarriers
    bestQAMSignal(num_sc,:) = antennaQAM(num_sc,:,bestAntenna(num_sc));
end


QAMRXSignalBest = reshape(bestQAMSignal,[],1); % reshape to QAM vector

% demod QAM
binaryRXDataBest = qamdemod(QAMRXSignalBest / sqrt(transmitterAveragePower), constellationCount,'OutputType','bit','UnitAveragePower',true);

% calculate error
error = sum(abs(binaryData - binaryRXDataBest),'all') / (bitPerConstellation * symbolsPerFrame * ( dataSubcarriers ));
fprintf("Best Antenna Sample, Bit Error: %.3f%%, SNR: %.2f [dB]\n", error*100, mean(mean(bestSNR)));

%% Stage 5: Display Metrics
% Display constellation map of first antenna signal (yellow) versus
% best-sampled signal (blue)
lim = sqrt(transmitterAveragePower)*3/2;

cd = comm.ConstellationDiagram(...
    'ShowReferenceConstellation',false,...
    'XLimits',[-lim lim],'YLimits',[-lim lim]...
);

step(cd,[ reshape(antennaQAM(:,:,1), [], 1) reshape(bestQAMSignal, [], 1) ]  )

% Display SNR of each antenna compared to best SNR from selective sampling
figure
hold on

for i = 1:antennaCount
    scatter(dataIndices, antennaSNR(:,i), 25, 'red', 'filled')
end
scatter(dataIndices, bestSNR, 15, 'green', 'filled')
hold off

grid on
title(sprintf("SNR Per OFDM Subcarrier Index Per Receiver (%d Antennas)", antennaCount))
ylabel("SNR [dB]")
xlabel("Subcarrier Index")

% Display SNR when the antenna count is increased
figure
hold on

antennaSampleCount = 1:antennaCount;
antennaSampleSNR = 1:antennaCount;
for i = antennaSampleCount
    antennaSampleSNR(i) = mean( max(antennaSNR(:,1:i),[],2) );
end
plot(antennaSampleCount, antennaSampleSNR,'*-','MarkerSize',10,'MarkerEdgeColor','red');

grid on
title(sprintf("Mean SNR per Number of Antenna Samples (Up to %d Antennas)", antennaCount))
ylabel("SNR [dB]")
xlabel("Number of Antenna Samples");
