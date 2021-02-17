% Alban-Félix Barreteau, M1 CORO SIP
% Mail : alban-felix.barreteau@eleves.ec-nantes.fr
% Matlab R2020a Update 5 (9.8.0.1451342), Student license
% Signal Filtering and System Identification (SISIF)
% Labwork 1 : Audio signal filtering

% Professor : said.moussaoui@ec-nantes.fr

%----------------------------------------------------------------%

clc

% Definition of your screen :
hpixels=1440;
vpixels=900;

%% Signal generation and analysis

% Signal simulation
duration=1; %duration of each note
fs=8192; %sampling frequency
fmin=250; fmax=550; f = linspace(fmin,fmax,1000);
[sig,t]=gamme_lab1(duration, fs); %simulate the signal
%soundsc(sig, fs); %play the signal

%Spectral analysis of the signal
fsig = [0:length(sig)-1]*(fs/length(sig));
S = fft(sig)/fs;

%Plotting the temporal signal and its spectrum
figure(10), hold off
set(figure(10),'name','Temporal gam signal and its spectrum','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])
subplot(2,1,1)
plot(t,sig)
ymin=-1.1*max(abs(sig)); ymax=1.1*max(abs(sig)); axis([0 8 ymin ymax])
xlabel('Frequency (in Hz)'), ylabel('Power')

subplot(2,1,2)
plot(fsig, abs(S)); 
xmin=fmin; xmax=fmax; ymin=0; ymax=1.1*max(abs(S)); axis([xmin xmax ymin ymax])
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal',"location",'southeast')



%% The retained specifications of the analog low-pass filter

fc=420; %The cutoff frequency
deltaf=20; %The transition band's width frequency
fpb=fc-deltaf/2; %Passband frequency
fsb=fc+deltaf/2; %Stopband frequency
Wp=2*pi*fpb; %Normalized passband edge frequency
Ws=2*pi*fsb; %Normalized stopband edge frequency
Rp=1; %Ripple in the passband
Rs=40; %Ripple in the stopband



%% Required order of the 4 filters

% Butterworth filter, N the lowest order and Wn the Butterworth natural frequency
[NButterworth, WnButterworth]=buttord(Wp,Ws,Rp,Rs,'s'); %for an analog filter, Wp and Ws in rad/s
[numButterworth,denButterworth]=butter(NButterworth,WnButterworth,'s');
HcButterworth=freqs(numButterworth,denButterworth,2*pi*f);
disp(['The required order for the Butterworth filter is ',num2str(NButterworth)])

% Chabyshev filter order 1, N the lowest order and Wp the Chebyshev natural frequency
[NChabyshev1, WpChabyshev1]=cheb1ord(Wp, Ws, Rp, Rs, 's'); %for an analog filter, Wp and Ws in rad/s
[numChabyshev1,denChabyshev1]=cheby1(NChabyshev1,Rp,WpChabyshev1,'s');
HcChabyshev1=freqs(numChabyshev1,denChabyshev1,2*pi*f);
disp(['The required order for the Chabyshev order 1 filter is ',num2str(NChabyshev1)])


% Chabyshev filter order 2, N the lowest order and Wp the Chebyshev natural frequency
[NChabyshev2, WpChabyshev2]=cheb2ord(Wp, Ws, Rp, Rs, 's'); %for an analog filter, Wp and Ws in rad/s
[numChabyshev2,denChabyshev2]=cheby2(NChabyshev2,Rs,WpChabyshev2,'s');
HcChabyshev2=freqs(numChabyshev2,denChabyshev2,2*pi*f);
disp(['The required order for the Chabyshev order 2 filter is ',num2str(NChabyshev2)])


% Cauer filter, N the lowest order and Wp the elliptic natural frequency
[NCauer, WpCauer] = ellipord(Wp, Ws, Rp, Rs, 's');
[numCauer, denCauer]=ellip(NCauer,Rp,Rs,Wp,'s');
HcCauer=freqs(numCauer,denCauer,2*pi*f);
disp(['The required order for the Cauer filter is ',num2str(NCauer)])




%% Magnitude frequency response of the 4 filters

figure(20), hold off
set(figure(20),'name','Magnitude frequency responses of the 4 filters','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])
plot(f,20*log10(abs(HcButterworth))); hold on
plot(f,20*log10(abs(HcChabyshev1)));
plot(f,20*log10(abs(HcChabyshev2)));
plot(f,20*log10(abs(HcCauer)));
xlabel('Frequency (in Hz)'), ylabel('Magnitude (in dB)'), grid on
legend('Butterworth','Chabyshev order 1','Chabyshev order 2','Cauer',"location",'southeast')
% The most advantageous filter is the cheaper one wich means the one with the lowest order



%% Signal filtering and spectral analysis

% Butterworth filter
sigfButterworth = lsim(tf(numButterworth,denButterworth),sig,t); %soundsc(sigfButterworth, fs);
fsigButterworth=[0:length(sigfButterworth)-1]*(fs/length(sigfButterworth));
SButterworth=fft(sigfButterworth)/fs;

figure(31), hold off
set(figure(31),'name','Original signal and Butterworth filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfButterworth);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Butterworth filtered signal',"location",'southeast')

subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigButterworth,abs(SButterworth))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Butterworth filtered signal',"location",'southeast')
try
axis([fmin fmax min(abs(SButterworth)) 1.1*max(abs(SButterworth))])
catch disp('Error during the numerical computation')
end


% Chabyshev filter of order 1
sigfChabyshev1 = lsim(tf(numChabyshev1,denChabyshev1),sig,t); %soundsc(sigfChabyshev1, fs);
fsigChabyshev1=[0:length(sigfChabyshev1)-1]*(fs/length(sigfChabyshev1));
SChabyshev1=fft(sigfChabyshev1)/fs;

figure(32), hold off
set(figure(32),'name','Original signal and Chabyshev1 filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfChabyshev1);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Chabyshev1 filtered signal',"location",'southeast')

subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigChabyshev1,abs(SChabyshev1))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Chabyshev1 filtered signal',"location",'southeast')
axis([fmin fmax min(abs(SChabyshev1)) 1.1*max(abs(SChabyshev1))])



% Chabyshev filter of order 2
sigfChabyshev2 = lsim(tf(numChabyshev2,denChabyshev2),sig,t); %soundsc(sigfChabyshev2, fs);
fsigChabyshev2=[0:length(sigfChabyshev2)-1]*(fs/length(sigfChabyshev2));
SChabyshev2=fft(sigfChabyshev2)/fs;

figure(33)
set(figure(33),'name','Original signal and Chabyshev2 filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfChabyshev2);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Chabyshev2 filtered signal',"location",'southeast')

subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigChabyshev2,abs(SChabyshev2))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Chabishev2 filtered signal',"location",'southeast')
axis([fmin fmax min(abs(SChabyshev2)) 1.1*max(abs(SChabyshev2))])



% Cauer filter
sigfCauer = lsim(tf(numCauer,denCauer),sig,t); %soundsc(sigfCauer, fs);
fsigCauer=[0:length(sigfCauer)-1]*(fs/length(sigfCauer));
SCauer=fft(sigfCauer)/fs;

figure(34)
set(figure(34),'name','Original signal and Cauer filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfCauer);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Cauer filtered signal',"location",'southeast')

subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigCauer,abs(SCauer))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Cauer filtered signal',"location",'southeast')
axis([fmin fmax min(abs(SCauer)) 1.1*max(abs(SCauer))])



%% The retained specification for the digital filter

fcl=340; %The low cutoff frequency
fch=360; %The high cutoff frequency
deltafl=10; %The low transition band's width frequency
deltafh=10; %The high transition band's width frequency
Rp=1; %Ripple in the passband
Rs=40; %Ripple in the stopband


%% Digital RII filter Chabyshef of order 2
% fpb and fsb should have 2 components because we design a stopband filter

fpl=fcl-deltafl/2; %Low passband edge frequency
fph=fch+deltafh/2; %High passband edge frequency
fpb=[fpl fph]; %Passband edges frequency
Wp=2*fpb/fs; %Normalized passband edge frequency

fsl=fcl+deltafl/2; %Low stopband edge frequency
fsh=fch-deltafh/2; %High stopband edge frequency
fsb=[fsl fsh]; %Stopband edges frequency
Ws=2*fsb/fs; %Normalized stopband edge frequency

% Chabyshev filter order 2, N the lowest order and Wp the Chebyshev natural frequency
[NChabyshev2, WpChabyshev2]=cheb2ord(Wp, Ws, Rp, Rs); %for a digital filter
[numChabyshev2,denChabyshev2]=cheby2(NChabyshev2,Rs,WpChabyshev2,'stop');
HdChabyshev2=freqz(numChabyshev2,denChabyshev2,f,fs);
disp(['The required order for the Chabyshev order 2 digital filter is ',num2str(NChabyshev2)])


figure (40), hold off 
set(figure(40),'name','Magnitude frequency response of the Chabyshev2 digital stopband filter','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])
plot(f,20*log10(abs(HdChabyshev2))) %Frequency response of the RII filter 
xlabel('Frequency (in Hz)'), ylabel('Magnitude (in dB)'), grid on
legend('Digital Chabyshev order 2 filter',"location",'southeast')




% Signal filtering and spectral analysis
sigfChabyshev2=filter(numChabyshev2,denChabyshev2,sig); %soundsc(sigfChabyshev2, fs);
fsigChabyshev2=[0:length(sigfChabyshev2)-1]*(fs/length(sigfChabyshev2));
SChabyshev2=fft(sigfChabyshev2)/fs;

figure(50)
set(figure(50),'name','Original signal and digital Chabyshev2 filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfChabyshev2);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Digital Chabyshev2 filtered signal',"location",'southeast')

subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigChabyshev2,abs(SChabyshev2))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Digital Chabishev2 filtered signal',"location",'southeast')
axis([fmin fmax min(abs(SChabyshev2)) 1.1*max(abs(SChabyshev2))])


%% FIR filter having the lowest lenght

% Calculation of the filter length
fcuts=[fpl fsl fsh fph]; %Specifications of characteristics frequencies sorted in increasing values
mags=[1 0 1]; %Filter magnitude in each band (linear scale)
epsilonp=1-10^(-Rp/20);
deltaa=10^(-Rs/20);
devs=[epsilonp deltaa epsilonp]; %Magnitude of admitted ripples in the bands for a stopband filter
[N,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fs);
disp(['The required order for the FIR digital filter is ',num2str(N)])


% Calculation of the FIR filter by the window method
window='kaiser';
numFIR=fir1(N,Wn,ftype,kaiser(N+1,beta)); %Window based FIR filter design
HdFIR=freqz(numFIR,1,f,fs);

figure(60), hold off
plot(f,20*log10(abs(HdFIR))); %Frequency response of the filter
xlabel('Frequency (in Hz)'), ylabel('Magnitude (in dB)'), grid on
legend('Digital FIR filter with the lowest lenght, Kaiser window',"location",'southeast')
set(figure(60),'name','Frequency response of the FIR filter','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])


% Signal filtering and spectral analysis
sigfFIR=filter(numFIR,1,sig); %soundsc(sigfFIR,fs)
fsigFIR=[0:length(sigfFIR)-1]*(fs/length(sigfFIR));
SFIR=fft(sigfFIR)/fs;

figure(70)
subplot(2,1,1)
plot(t,sig); hold on; plot(t,sigfFIR);
xlabel('Time (in s)'), ylabel('Magnitude')
legend('Original signal','Digital FIR filtered signal',"location",'southeast')
subplot(2,1,2)
plot(fsig, abs(S)), hold on; plot(fsigFIR,abs(SFIR))
xlabel('Frequency (in Hz)'), ylabel('Power')
legend('Original signal','Digital FIR filtered signal',"location",'southeast')
axis([fmin fmax min(abs(SFIR)) 1.1*max(abs(SFIR))])
set(figure(70),'name','Original signal and digital FIR filtered signal','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])


%% FIR and RII filters comparison

figure(80), hold off
plot(f,20*log10(abs(HdFIR))), hold on %Frequency response of the FIR filter
plot(f,20*log10(abs(HdChabyshev2))) %Frequency response of the RII filter
legend('FIR filter with the lower lenght, Kaiser window','RII Chabyshev order 2 filter','location','southeast'), grid on
xlabel('Frequency (in Hz)'); ylabel('Magnitude (in dB)');
set(figure(80),'name','FIR and RII filter magnitude frequency responses comparison','position',[hpixels/2 vpixels/2 hpixels/2 vpixels/2])

