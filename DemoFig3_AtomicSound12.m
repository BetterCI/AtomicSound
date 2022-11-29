clear 


% 2022/11/29 
% Author: F. Kong, Q. Meng. Email: mengqinglin@scut.edu.cn
% Use this program to generate Fig.3 of the manuscript of
% "A generic signal processing framework for speech redundancy manipulation algorithms in speech perception studies"
% 

%% Initial parameters
soundRootPath = ''; corpusPath = 'sounds\'; filename = 'mhint001.wav';
soundPath = [soundRootPath,corpusPath,filename];

[x0,fs] = audioread(soundPath);

if fs~=16000
    x0 = resample(x0,16000,fs); fs = 16000;
end
x0 = x0 * 0.3;
figure;
t = (0:length(x0)-1)/fs;
subplot(621);plot(t,x0);set(gca,'XTick',[])
subplot(6,2,[3,5]);myspectrogram(x0,fs);set(gca,'XTick',[]);xlabel('')

for demo_vocoder_type = [1,2,3] % number of maxima  1. Nmax = 1; 2. Nmax = 4; 3. Nmax = 8
    clear p
    switch demo_vocoder_type
        case 1
            p.num_Ch = 32;  % channel number,m
            p.max_Ch_num = 1;
            p.rat_pCh = 106; % Pulse rate, pps
      
        case 2
            p.num_Ch = 32;  % channel number,m
            p.max_Ch_num = 4;
            p.rat_pCh = 18; % Pulse rate, pps                        
        case 3
            p.num_Ch = 32;  % channel number,m
            p.max_Ch_num = 8;
            p.rat_pCh = 13; % Pulse rate, pps
        
    end

            p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking
            p.carr_typ = 1; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                            % 4. TFS from another sound,i.e., chameric sounds  

    p.pre_emp = 1;
    p.bfilt_typ = 2; % gammatone filter bank  只适用于30个以上的通道数; 2 butterworth filter bank
    p.fre_range = [80, 7990];% Hz
    p.env_cut = 200;% Hz
        if p.max_mod >1 && p.max_mod <5
            p.block_T = 0.1;% seconds (useful only when p.max_mod > 1)
        elseif p.max_mod == 5
            p.bubbleNum = 8;
            p.maskingFlag = 0;% 0 区域内删除  1 区域内为保留(暂时不可用，还没调试好)；
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
    switch demo_vocoder_type
        case 1
            subplot(6,2,2);plot(t,Output1);set(gca,'XTick',[])
            subplot(6,2,[4,6]);myspectrogram(Output1,fs);set(gca,'XTick',[]);xlabel('')
        case 2
            subplot(6,2,7);plot(t,Output1);set(gca,'XTick',[])
            subplot(6,2,[9,11]);myspectrogram(Output1,fs);
        case 3
            subplot(6,2,8);plot(t,Output1);set(gca,'XTick',[])
            subplot(6,2,[10,12]);myspectrogram(Output1,fs);
    end

    audiowrite(['demo_atomic_',num2str(demo_vocoder_type),'.wav'],Output1,fs);
    %     soundsc(Output,fs);
end
    a = findobj('Type','Axes');linkaxes(a,'x');

