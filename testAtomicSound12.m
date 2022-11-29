clear 

%% Initial parameters
p.pre_emp = 1;
p.num_Ch = 32;  % channel number,m
p.bfilt_typ = 2; % gammatone filter bank  for num_Ch>30; 2 butterworth filter bank
p.max_Ch_num = 2;
p.fre_range = [80, 7990];% Hz
p.env_cut = 200;% Hz
p.rat_pCh = 100; % Pulse rate, pps
p.max_mod = 1; % 1 n-of-m; 2 net_mosaic 3 env rms equalization mosaic 4. Time reversal Envelope 5. bubble masking

p.tim_typ = 4; 
% 1 randomized within each frame for each band; 2 low to high band; 3 high to low band ; 4 Peak time in each channel-frame
% 
p.jit_max = 0;
p.carr_typ = 1; % 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise;
                % 4. TFS from another sound,i.e., chameric sounds
%p.spread = -12;  % dB/oct default is -12. if you don't want to add interchannel spread, please comment it out.   
% p.quantization = {8,45}; % bits and dynamic range

    if p.max_mod >1 && p.max_mod <5
        p.block_T = 0.05;% seconds (useful only when p.max_mod > 1)
    elseif p.max_mod == 5
        p.bubbleNum = 8;
        p.maskingFlag = 0;%
        p.bubble_D = [0.1,0.3];  % the bubble duration will be a value which is uniformly distributed in this range. if you want a constant duration, you can make the two values the same.
        p.bubble_B = [0.4,0.8];  % the bubble bandwidth will be a value which is uniformaly distributed the range of center frequency * p.bubble_B
    end
    
%% load wav
% soundRootPath = 'D:\นคื๗\Library\Sound\'; corpusPath = 'MSP_FU\normal16kHz\'; filename = 'msp019.wav';
soundRootPath = ''; corpusPath = 'sounds\'; filename = '203.wav';
soundPath = [soundRootPath,corpusPath,filename];


[x0,fs] = audioread(soundPath);
if length(x0)/fs > 3
    x0 = x0(1:3*fs);
end

if fs~=16000
    x0 = resample(x0,16000,fs); fs = 16000;
end

if p.carr_typ == 4 % for chameric
    sound2filename = 'Adult004';
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



Output = [];
for n = 1:size(x0,2)
    x = x0(:,n);  
    [Output1,GEV_sound_Ch] = AtomicSound12(x,fs,p);

    t = (0:length(Output1)-1)/fs;
    GEV_sound_Ch = GEV_sound_Ch / max(GEV_sound_Ch(:)) * max(Output1);
    figure(2*n-1);
    plot(t,Output1);

    %%
%     T = [5,5];
    figure(2*n);
    subplot(711);plot(t,x(1:length(t))); axis tight; 

    subplot(7,1,[2,3]);myspectrogram(x,fs);
    subplot(714); 
    plot(t,Output1);
    subplot(715); 
    
    for k = 1:p.num_Ch
        tempChSig = GEV_sound_Ch(k,1:length(t));
        abs_tempChSig = abs(tempChSig);
        tempChSig(abs_tempChSig<1e-5) = NaN;
         plot(t,tempChSig,'Color',[k/p.num_Ch-0.4*k/p.num_Ch,1-0.9*k/p.num_Ch,0]); hold on;
    end
    axis tight; hold off
    subplot(7,1,[6,7]);myspectrogram(Output1,fs);
    Output = [Output, Output1];
end
    soundsc([Output],fs);

a = findobj('Type','Axes');linkaxes(a,'x');

audiowrite('Output.wav',Output,fs)