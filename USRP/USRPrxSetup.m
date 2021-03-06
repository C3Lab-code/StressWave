decimationfactor = 100e6/fs;
rx = comm.SDRuReceiver(...
              'Platform','N200/N210/USRP2', ...
              'IPAddress','192.168.10.2', ...
              'CenterFrequency',0, ...
              'EnableBurstMode', true,...
              'SamplesPerFrame', 6000*2,...
              'NumFramesInBurst', 10,...
              'DecimationFactor',decimationfactor);