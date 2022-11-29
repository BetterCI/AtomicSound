# AtomicSound
MATLAB program for a generic framework of speech degradation algorithms and its special case atomic speech


To make the framework and atomic speech model available to public for research purpose, the program is provided here freely for acadamic research.
Citation:


Fanhui Kong, Huali Zhou, Qinglin Meng, Nengheng Zheng. "A generic signal processing framework for speech redundancy manipulation algorithms in speech perception studies"
Submitted to a journal.


Main program: AtomicSound12.m
Demo1: testAtomicSound12.m
Demo2: DemoFig2_AtomicSound12.m  for Fig.2 in the manuscript
Demo3: DemoFig3_AtomicSound12.m  for Fig.3 in the manuscript

Funtion: a signal processing framework that can generate several degradation algorithms by simply switching several key parameters
         The algoirthms include mosaic speech, chimaeric speech, sine-wave carrier and noise carrier channel vocoders for cochlear implant simulations, and pulsatile Gaussian-enveloped-Tones vocoder for cochlear implant simulation.
         
testAtomicSound12.m Version 12.0; Date: 2022/11/15; Email: mengqinglin@scut.edu.cn
Atomic sound is defined as the sound signal resynthesized using tiny sound atoms, i.e. Gaussian enveloped pulses
The original idea is from the Gaussian-Enveloped-tone vocoder which was proposed to simulate the pulsatile stimulation of cochlear implants. 
The first ever proposal for using the Gaussian-Enveloped-Tone to study hearing may be Gabor(1942) which borrow ideas from the quantum theory. 
The uncertainty theory in both quantum theory and signal processing share the same mathmatical basis. 
The following is the parameter introduction
Inputs:
           in_sig          input signal
           fs              sampling rate (Hz)
           p               paramters
           p.pre_emp       0. Without pre-emphasis; 1(default). With pre-emphasis
           p.num_Ch        number of frequency channels
           p.fre_range     analysis frequency range (Hz)
           p.env_cut       the cutoff frequency for envelope extraction within each channel (Hz)
           p.rat_pCh       pulse rate per channel (pps or Hz)
           p.max_Ch_num    number of selection of the channels with maximum information
           p.max_mod       mode of selection of the channels with maximum information. The default criterion is the channels with maximum intensity
           p.bfilt_typ     type of bandpass filter bank, 
                           1. Gammatone, only used for >30 num_Ch, 2. Butterworth
           p.tim_typ       the type of sampling point
                                   1 randomized within each frame for each band; 
                                   2 low to high band; 
                                   3 high to low band; 
                                   4 Peak in each channel-frame
           p.jit_max       maximum jitter time (s); default = 0 s; may be used when p.tim_typ = 2 or 3 to reduce the regularity? 
           p.carr_typ      carrier type; 1.Origial temporal fine structure (TFS) from the in_sig; 2. sine wave; 3. band-limited noise
           p.block_T(optional)      the duration of mosaic blocks
           p.Spread        dB/oct