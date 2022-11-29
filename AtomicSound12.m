function [out_sig,GEV_sound_Ch] = AtomicSound12(in_sig,fs,p)
% The function is used to analysis the input sound and resynthesize it into a pulsatile (or atomic) sound
% Author: Fanhui Kong, Qinglin Meng
% Version 9.0; Date: 2022/4/20; Email: mengqinglin@scut.edu.cn
% Citation: Meng, Q. (2020). Perception of atomic speech. The Journal of the Acoustical Society of America, 148(4), 2722-2722.

% To serve the framework study
% Version 10.0; Date: 2022/9/12; Email: mengqinglin@scut.edu.cn

% add more 
% Version 11.0; Date: 2022/10/21; Email: mengqinglin@scut.edu.cn

% remove spectral smearing from the step before temporal sparsing or smearing to the step before the modulation (AM/FM)
% Version 12.0; Date: 2022/11/15; Email: mengqinglin@scut.edu.cn
% Citation: F. Kong, H. Zhou, N. Zheng, and Q. Meng. (2020). 
% A generic signal processing framework for speech reduandacy manipulation algorithms in speech perception studies.

% Atomic sound is defined as the sound signal resynthesized using tiny sound atoms, i.e. 
% Gaussian enveloped pulses
% The original idea is from the Gaussian-Enveloped-tone vocoder which was
% proposed to simulate the pulsatile stimulation of cochlear implants. 
% The first ever proposal for using the Gaussian-Enveloped-Tone to study
% hearing may be Gabor(1942) which borrow ideas from the quantum theory. 
% The uncertainty theory in both quantum theory and signal processing share
% the same mathmatical basis. 
% The following is the parameter introduction
% % Inputs:
%           in_sig          input signal
%           fs              sampling rate (Hz)
%           p               paramters
%           p.pre_emp       0. Without pre-emphasis; 1(default). With pre-emphasis
%           p.num_Ch        number of frequency channels
%           p.fre_range     analysis frequency range (Hz)
%           p.env_cut       the cutoff frequency for envelope extraction within each channel (Hz)
%           p.rat_pCh       pulse rate per channel (pps or Hz)
%           p.max_Ch_num    number of selection of the channels with maximum information
%           p.max_mod       mode of selection of the channels with maximum
%                           information. The default criterion is the channels 
%                           with maximum intensity
%           p.bfilt_typ     type of bandpass filter bank, 
%                           1. Gammatone, only used for >30 num_Ch, 2. Butterworth
%           p.tim_typ       the type of sampling point
%                                   1 randomized within each frame for each band; 
%                                   2 low to high band; 
%                                   3 high to low band; 
%                                   4 Peak in each channel-frame
%           p.jit_max       maximum jitter time (s); default = 0 s; may be used when p.tim_typ = 2 or 3 to reduce the regularity? 
%           p.carr_typ      carrier type; 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise
%           p.block_T(optional)      the duration of mosaic blocks
%           p.Spread        dB/oct
% % Outputs:

targetRMS = rms(in_sig);
%% pre-empphasis
if ~isfield(p,'pre_emp') || p.pre_emp
    [bb,aa]=butter(1,1200/(fs/2),'high'); 
    x=filter(bb,aa,in_sig);
    x = x * targetRMS / rms(x);
else
    x = in_sig;
end
t = (0:length(x)-1)/fs;
% SR_Pulses = p.SR_Pulses; % Pulse rate, pps
% numChan = p.numChan;  % channel number,m
% ZeroCh = p.ZeroCh;
% synthesisFilterBankType = p.synthesisFilterBankType;% 1 Gammatone 一般适用于30个通道以上, 2 Butterworth
% tfsOrNot = p.TFSFlag; % 1 包含TFS，不包含TFS
% frequencyMappingShift = p.frequencyMappingShift;
% fRange = p.fRange;
% pulseTimingFlag = p.pulseTimingFlag; % 1 randomized within each frame for each band; 2 low to high band; 3 high to low band ; 4 Peak in each channel-frame
% maxTimeJitter = 0;
%% Intial parameter check
if p.num_Ch < 30 && p.bfilt_typ == 1
    disp('Warning: channel number < 30 is not suggetted to be used with Gammatone filter bank')
end

%% Step 1： Band-pass filter bank
[cf,fcs,r,b,a] = mybandpassfilterbank(x,fs,p.num_Ch,p.fre_range,p.bfilt_typ);

%% Step 2: Envelope extraction
analyticSignal = (hilbert(r'))'; h_env_perio = abs(analyticSignal); % envelope and periodicity

env_cut = min([p.env_cut,p.rat_pCh/2]);
[b_env,a_env] = butter(3,env_cut/(fs/2)); % 提取包络，用于进行n-of-m最大值选取
for n = 1:p.num_Ch  % 
    h_env(n,:) = filtfilt(b_env,a_env,h_env_perio(n,:)); % smoothing the envelope
end

perio_cut = p.rat_pCh/2;
[b_perio,a_perio] = butter(3,perio_cut/(fs/2));  % 按照刺激速率保留一定的周期性性信息。用于与TFS合成。
for n = 1:p.num_Ch  % 
    h_env_perio(n,:) = filtfilt(b_perio,a_perio,h_env_perio(n,:)); % smoothing the envelope
end

h_env_perio(h_env_perio<0) = 1e-10;% 20221029 如果非正则置为非常小的正值


%% Step 5: Intensity sparsing
if isfield(p,'quantization')
    sort_h_env_perio = sort(h_env_perio(:));
    ref_h_env_perio = sort_h_env_perio(round(length(sort_h_env_perio)*0.95)); % 从大到小排序，其中排在95%的较大值作为参考
    h_env_perio(h_env_perio>ref_h_env_perio) = ref_h_env_perio; %wave clipping
    dB_h_env_perio = 20*log10(h_env_perio/ref_h_env_perio);
    dynamicrange = cell2mat(p.quantization(2));
    targetbits = cell2mat(p.quantization(1));
    dB_h_env_perio = dB_h_env_perio + dynamicrange;
    dB_h_env_perio(dB_h_env_perio<0) = 0;
    quantized_h_env_perio = fix(dB_h_env_perio / dynamicrange * (2^targetbits-1))/((2^targetbits-1))*dynamicrange-dynamicrange;
    h_env_perio = 10.^(quantized_h_env_perio/20)*ref_h_env_perio;
%     h_env_perio(dB_h_env_perio<0) = 0;
end



%% Step 2 TFS Extraction
% TFS 
h_tfs = zeros(size(h_env));
temp_tfs = h_tfs;
switch p.carr_typ 
    case 1 % original tfs
        h_tfs = cos(angle(analyticSignal));
        for n = 1:p.num_Ch
           temp_tfs(n,:) = filtfilt(b(n,:),a(n,:),h_tfs(n,:));
           h_tfs(n,:) = temp_tfs(n,:) / rms(temp_tfs(n,:)) * rms(h_tfs(n,:));
        end
    case 2 % sine wave
        for n = 1:p.num_Ch
            h_tfs(n,:) = sin(2*pi*(cf(n))*t+rand(1)*2*pi);
        end
    case 3 % noise carrier
        for n = 1:p.num_Ch
%             bandwidth = diff(cf);
            sin_component_num = ceil(cf(n)*0.1);
            for m = 1:sin_component_num % 
                tempFreq = rand(1)*(fcs(n+1)-fcs(n))+fcs(n);
                h_tfs(n,:) = h_tfs(n,:)+sin(2*pi*tempFreq*t+rand(1)*2*pi);% adding many sines with random initial phase gets a noise
            end
            h_tfs(n,:) = h_tfs(n,:) / sin_component_num;
        end
    case 4 % chimeric sound
        %% Band-pass filter bank, env/tfs decomposition
        if length(p.sound2) ~= length(x)
                disp('the two sounds for chameric synthesis have different length');
        end
        [~,~,r_sound2,b,a] = mybandpassfilterbank(p.sound2,fs,p.num_Ch,p.fre_range,p.bfilt_typ);
        analyticSignal_sound2 = (hilbert(r_sound2'))'; %这两句应该与69 70 行对in_sig的处理保持一致
        
        h_tfs = cos(angle(analyticSignal_sound2));
        for n = 1:p.num_Ch
           temp_tfs(n,:) = filtfilt(b(n,:),a(n,:),h_tfs(n,:));
           h_tfs(n,:) = temp_tfs(n,:) / rms(temp_tfs(n,:)) * rms(h_tfs(n,:));
        end        
end

%% Step 3 and 4: Temporal and Spectral sparsing or smearing

framelength = floor(fs/p.rat_pCh);
frameNum = floor(length(x)/framelength);
newh = zeros(size(h_env));
for n = 1:p.num_Ch
    switch p.tim_typ
        case 1
            nChanPoints = (1:framelength:framelength*frameNum)+randi(framelength,1,frameNum)-1;
        case 2
            nChanPoints = (1:framelength:framelength*frameNum)+floor(framelength/p.num_Ch*(n-1));
        case 3
            nChanPoints = (1:framelength:framelength*frameNum)+floor(framelength-framelength/p.num_Ch*(n-1));
        case 4
            for m = 1:frameNum
                tempindex = ((m-1)*framelength+1):m*framelength;
                tempPeakindex = (m-1)*framelength + find(h_env(n,tempindex) == max(h_env(n,tempindex)));
                if length(tempPeakindex) > 1
                    tempPeakorder = randperm(length(tempPeakindex));
                    tempPeakindex = tempPeakindex(tempPeakorder(1));
                end
                nChanPoints(m) = tempPeakindex;
            end
    end
    ChanPoints(n,:) = nChanPoints;
    egram(n,1:frameNum) = h_env_perio(n,nChanPoints);
end



newegram = egram;
switch p.max_mod
    case 1 % channels with maxima intensities, i.e., n-of-m used in some cochlear implant strategies.
        if p.max_Ch_num < p.num_Ch
            ZeroCh = 1:p.num_Ch-p.max_Ch_num;
        else
            ZeroCh = [];
        end
        for m = 1:frameNum
            [a1,b1] = sort(newegram(:,m));
             zeroChs = b1(ZeroCh);
             newegram(zeroChs,m) = 0;
        end     
    case 2 % 按照网格状进行切分  默认采用fcs为横边界，输入参数block_T的整数倍为纵边界
            
        for m = 1:frameNum
            mosaicBlock = ceil(m*framelength/fs/p.block_T);
            if mod(mosaicBlock,2) % 奇数
                zeroChs = 1:2:p.num_Ch;
            else % 偶数
                zeroChs = 2:2:p.num_Ch;
            end
            newegram(zeroChs,m) = 0;
        end
    case 3 % 按照网格状进行切分  默认采用fcs为横边界，输入参数block_T的整数倍为纵边界,并且进一步把网格内的包络平均

        mosaicBlocks = ceil((1:frameNum)*framelength/fs/p.block_T);
        mosaicBlocksNum = mosaicBlocks(end);
        for m = 1:mosaicBlocksNum
            for n = 1:p.num_Ch
                temp_index = mosaicBlocks==m;
                averg_env_temp = rms(newegram(n,temp_index));
                newegram(n,temp_index) = ones(1,length(newegram(n,temp_index))) * averg_env_temp;
            end
        end       
    case 4 % 按照网格状进行切分  默认采用fcs为横边界，输入参数block_T的整数倍为纵边界,并且进一步把网格内的包络在时间上反转
        mosaicBlocks = ceil((1:frameNum)*framelength/fs/p.block_T);
        mosaicBlocksNum = mosaicBlocks(end);
        for m = 1:mosaicBlocksNum
            temp_index = find(mosaicBlocks==m);
            newegram(:,temp_index) = newegram(:,fliplr(temp_index));
        end     
    case 5 % bubble masking, 生成若干个圆形掩蔽区域或黑白颠倒
        for m = 1:p.bubbleNum 
            % each a bubble is assummed to be an ellipse whose stardard equation is x^2/a^2 + y^2/b^2 = 1
            bubbleOrigin_T(m) = rand(1)*(frameNum-1)*framelength/fs;
            bubbleOrigin_F(m) = rand(1)*diff(p.fre_range) + p.fre_range(1);
            
            bubble_a(m) = rand(1)*diff(p.bubble_D) + p.bubble_D(1);
            bubble_b(m) = bubbleOrigin_F(m) * (rand(1)*diff(p.bubble_B) + p.bubble_B(1));
        end
        for n = 1:p.num_Ch
            for k = 1:frameNum
                temp_t_f = [(k-1)*framelength/fs,cf(n)];
                for m = 1:p.bubbleNum
                    if (temp_t_f(1)-bubbleOrigin_T(m))^2/(bubble_a(m)/2)^2 + (temp_t_f(2)-bubbleOrigin_F(m))^2/(bubble_b(m)/2)^2 < 1
                        newegram(n,k) = p.maskingFlag*newegram(n,k);
                    else 
                        newegram(n,k) = (1-p.maskingFlag)*newegram(n,k);
                    end
                end
            end
        end
end


for n = 1:p.num_Ch
    newh(n,ChanPoints(n,:)) = newegram(n,1:frameNum);
end


%% 
% cf = cf + frequencyMappingShift;

%%

if p.jit_max 
    maxTimeJitterSample = p.jit_max * fs;
    newnewh = zeros(size(newh));
    for n = 1:numChan
        for m = 1:size(newh,2)
            if newh(n,m)
                tempJitterSample = ceil((rand(1)*2-1)*maxTimeJitterSample);
                if m+tempJitterSample>size(newh,2)
                    newM = size(newh,2);
                elseif m+tempJitterSample<1
                    newM = 1;
                else
                    newM = m+tempJitterSample;
                    newnewh(n,newM) = newh(n,m) + newnewh(n,newM);
                end
            end
        end
    end

    newh = newnewh;
end
%% Step 5: Convoluation with Gaussian pulses

b = 1.019*24.7*(4.37*cf(1)/1000+1);       % rate of decay or bandwidth
D(1) = 1/b;
halfEnvelopePointNum = round(D(1)*fs);

t = (0:size(newh,2)-1)/fs;
for n = 1:p.num_Ch
    b = 1.019*24.7*(4.37*cf(n)/1000+1);       % rate of decay or bandwidth
    D(n) = 1/b; % D is the effective duration of Gaussian Envelope
    GauEnvFir = GaussianEnvelope(D(n),halfEnvelopePointNum,fs);
    temp = conv(newh(n,:),GauEnvFir);
    tempHEP = round(halfEnvelopePointNum/2);
    GEV_AM_Ch(n,:) = temp(tempHEP:length(t)+tempHEP-1);
end


%% Step 7 Spectral smearing
% current spread among channels
if isfield(p,'spread')
    cfMatrix = repmat(cf',1,p.num_Ch); % 每个通道的中心频率，复制N份
    for n = 1:p.num_Ch
        SpreadMatrixdB(:,n) = abs(log2(cfMatrix(:,n)/cf(n)))*p.spread;% spread的默认值可能为-8dB/oct
    end

    SpreadMatrix = 10.^(SpreadMatrixdB / 20);
    GEV_AM_Ch = (sqrt((GEV_AM_Ch').^2 * SpreadMatrix.^2))';   % 参考 Oxenham 2014
end


%% Step 8: Modulation (AM/FM)
for n = 1:p.num_Ch
    GEV_sound_Ch(n,:) = GEV_AM_Ch(n,:).*h_tfs(n,:);
    if rms(GEV_sound_Ch(n,:)) ~= 0
        GEV_sound_Ch(n,:) = GEV_sound_Ch(n,:) * rms(r(n,:)) / rms(GEV_sound_Ch(n,:));
    end
end




%% Step 9: Sum output channels


Output = sum(GEV_sound_Ch);
out_sig = Output(:);



out_sig = out_sig * targetRMS / rms(out_sig(:));


function GauEnv = GaussianEnvelope(D,halfEnvelopePointNum,fs)
t = (-halfEnvelopePointNum:halfEnvelopePointNum)/fs;
GauEnv = exp(-pi*t.^2/D^2);








