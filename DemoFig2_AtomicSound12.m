clear 

% 2022/11/29 
% Author: F. Kong, Q. Meng. Email: mengqinglin@scut.edu.cn
% Use this program to generate Fig.2 of the manuscript of
% "A generic signal processing framework for speech redundancy manipulation algorithms in speech perception studies (L)"
% 

%% Initial parameters
soundRootPath = ''; corpusPath = 'sounds\'; filename = 'mhint001.wav';
soundPath = [soundRootPath,corpusPath,filename];

[x0,fs] = audioread(soundPath);

if fs~=16000
    x0 = resample(x0,16000,fs); fs = 16000;
end
figure;
subplot(231);myspectrogram(x0,fs);
for demo_vocoder_type = [1,2,4,5,3] % 1. Sine-wave Channel Vocoder 2. Noise Channel Vocoder
                       % 3. Mosaic Speech; 4. Chimaric speech; 5. Gaussian-enveloped tones
    clear p
    switch demo_vocoder_type
        case 1
            p.num_Ch = 8;  % channel number,m
            p.max_Ch_num = 8;
            p.rat_pCh = 1000; % Pulse rate, pps
            p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 2; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds        
        case 2
            p.num_Ch = 8;  % channel number,m
            p.max_Ch_num = 8;
            p.rat_pCh = 1000; % Pulse rate, pps
            p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 3; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds                          
        case 3
            p.num_Ch = 8;  % channel number,m
            p.max_Ch_num = 8;
            p.rat_pCh = 1000; % Pulse rate, pps
            p.max_mod = 3; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 3; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds         
        case 4
            p.num_Ch = 16;  % channel number,m
            p.max_Ch_num = 16;
            p.rat_pCh = 2000; % Pulse rate, pps
            p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 4; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds   
        case 5
            p.num_Ch = 22;  % channel number,m
            p.max_Ch_num = 8;
            p.rat_pCh = 900; % Pulse rate, pps
            p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 2; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds         
            p.quantization = {8,55}; % bits and dynamic range
    end



    p.pre_emp = 1;
    p.bfilt_typ = 2; % gammatone filter bank  ????????30??????????????; 2 butterworth filter bank
    p.fre_range = [80, 7990];% Hz
    p.env_cut = 200;% Hz
        if p.max_mod >1 && p.max_mod <5
            p.block_T = 0.1;% seconds (useful only when p.max_mod > 1)
        elseif p.max_mod == 5
            p.bubbleNum = 8;
            p.maskingFlag = 0;% 0 ??????????  1 ????????????(??????????????????????)??
            p.bubble_D = [0.1,0.3];  % the bubble duration will be a value which is uniformly distributed in this range. if you want a constant duration, you can make the two values the same.
            p.bubble_B = [0.4,0.8];  % the bubble bandwidth will be a value which is uniformaly distributed the range of center frequency * p.bubble_B
        end

    p.tim_typ = 4; 
    % 1 randomized within each frame for each band; 2 low to high band; 3 high to low band ; 4 Peak time in each channel-frame
    % 
    p.jit_max = 0;

    % p.spread = -12;  % dB/oct default is -12. if you don't want to add interchannel spread, please comment it out.   
    % p.quantization = {8,45}; % bits and dynamic range



    if p.carr_typ == 4 % for chameric
        sound2filename = 'mhint002.wav';
        sound2Path = [soundRootPath,corpusPath,sound2filename];

        [sound2,fs] = audioread(sound2Path);
        if fs~=16000
            sound2 = resample(sound2,16000,fs); fs = 16000;
        end
        if length(sound2) > length(x0)
            sound2 = sound2(1:length(x0));
        else
            sound2 = [sound2;zeros(length(x0)-length(sound2),1)];
        end
        p.sound2 = sound2;
    end




        [Output1,GEV_sound_Ch] = AtomicSound12(x0,fs,p);

        t = (0:length(Output1)-1)/fs;
    subplot(2,3,demo_vocoder_type+1);myspectrogram(Output1,fs);
    audiowrite(['demo_vocoder_type',num2str(demo_vocoder_type),'.wav'],Output1,fs);
    %     soundsc(Output,fs);
end
    a = findobj('Type','Axes');linkaxes(a,'x');

