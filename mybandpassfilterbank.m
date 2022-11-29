function [cf,fcs,r,b,a] = mybandpassfilterbank(x,fs,chNum,frange,ftype)
% 20191230
% Bandpass filter-bank
% Input: x      input sound, N×1
%        fs     sampling rate (Hz)
%        chNum  band number
%        frange frequency range (Hz)
% Ouput: cf     center frequencis (Hz), 1× (chNum )
%        fcs    cutoff frequencies (Hz), 1 × (chNum +1)
%        r      band signals, chNum × N
%        b,a    filter parameters
switch ftype
    case 1  
        [cf,r,b] = gammatone(x, chNum, frange, fs);
        a = ones(chNum,1);
        fcs(1) = frange(1);
        fcs(chNum+1) = frange(end);
        for n = 2:chNum
            fcs(n) = mean(cf(n-1:n));
        end
    case 2
        erb_b = hz2erb(frange);       % upper and lower bound of ERB
        erb = [erb_b(1):diff(erb_b)/(chNum):erb_b(2)];     % ERB segment
        fcs = erb2hz(erb);
        for n = 1:chNum
            [b(n,:),a(n,:)] = butter(3,fcs(n:n+1)/(fs/2));
            r(n,:) = (filtfilt(b(n,:),a(n,:),x))';
            cf(n) = mean(fcs(n:n+1));
        end
end

