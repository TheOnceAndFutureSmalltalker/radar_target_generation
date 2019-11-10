clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
v = -20  % constant velocity 20m/s towards ego car
R = 110  % initial distance of car from ego car

 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

% calculae bandwidth or B seep
c = 3*10^8;
d_res = 1
B = c / (2 * d_res); % bandwidthm or B sweep

% calculate chirp time
Rmax = 200 % max range is 200 m
Tchirp = 5.5 * 2 * Rmax / c; % chirp time

% calculate slope
slope = B / Tchirp

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                     
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
 
for i=1:length(t)         
  
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v * t(i);    % distance to target for current timestep
    td(i) = 2 * r_t(i) / c;   % time of travel for signal for current timestep
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2 * pi * (fc * t(i)+ (slope * t(i) * t(i) /2) ));
    Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + (slope * (t(i) - td(i)) * (t(i) - td(i)) / 2 ) ));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix, Nr);
Mix_fft = Mix_fft./Nr;

% Take the absolute value of FFT output
Mix_fft = abs(Mix_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix_fft = Mix_fft(1:Nr/2);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
plot(Mix_fft);
 
axis ([0 200 0 1]);


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Select the number of Training Cells in both the dimensions.
Tr = 6;
Td = 6;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5;
Gd = 2;

% offset the threshold by SNR value in dB
offset = 9;

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

% iterating across range, start at cell Gr+Tr+1 so that window fits
rstart = Gr+Tr+1;
% using only one side of the range dimension, Nr/2 from above,
% stop with enough cells for last window to fully fit
rstop = Nr/2-Gr-Tr;
% same for iterating doppler dimension, full window must fit
dstart = Td+Gd+1;
dstop = Nd-Gd-Td;

% calculate number of training cells in window as the full window 
% of cells minus the center gurad cells and CUT
num_training_cells = (2*Td+2*Gd+1)*(2*Tr+2*Gr+1) - (2*Gd+1)*(2*Gr+1);

% start sliding window across, then down
for r = rstart:rstop
    for d = dstart:dstop
        
        % noise accumulator variable for window
        noise_level = zeros(1,1);
        
        % sum signal over just the training portions of window
        % thus ignoring the guard cells and CUT
        for w_r = r-Tr-Gr:r+Tr+Gr
            for w_d = d-Td-Gd:d+Td+Gd
                if( w_r<(r-Gr) | w_r>(r+Gr) | w_d<(d-Gd) | w_d>(d+Gd) )
                    % sum linear representation of signal
                    noise_level = noise_level + db2pow(RDM(w_r,w_d));
                end
            end
        end
        
        % average the noise_level by dividing the sum by
        % the number of training cells in window
        avg_noise_level = noise_level/num_training_cells;
        
        % convert back to logarithmic value
        avg_noise_level_db = pow2db(avg_noise_level);
        
        % add offset to get threshold value
        threshold = avg_noise_level_db + offset;
        
        % if signal above threshold, reassign to 1
        % else, reassign to 0
        CUT = RDM(r,d);
        if(CUT > threshold)
            RDM(r,d) = 1;
        else
            RDM(r,d) = 0 ;
        end
        
    end
end



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
% set perimeter, unthresholded values, of RDM to zeroes
for r = 1:Nr/2
    for d = 1:Nd
        if( r<(Tr+Gr+1) | r>(Nr/2-Tr-Gr-1) )
            RDM(r,d) = 0;
        end
        if( d<(Td+Gd+1) | d>(Nd-Td-Gd-1) )
            RDM(r,d) = 0;
        end
    end
end



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM);
colorbar;


 
 