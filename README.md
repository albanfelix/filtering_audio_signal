# Audio signal filtering
This script gives an example of filter synthesis and its application to an audio signal processing (elimination or extraction of notes) using the available functions in Matlab.

# Filtering of an audio signal
Musical notes of a piano can be simulated numerically by assuming that each note corresponds to a pure sine signal with frequency value are given in table. A Matlab function gamme.m allowing to simulate a sequence of eight notes successively is given.
Note Do Re Mi Fa Sol La Si Do (Frequency(Hz) 262   294  330   349   392   440   494   523)
Signal generation and analysis : Generate the signal and listen it. The duration is one second for each note and the sampling frequency is Fs = 8192 Hz. Plot the temporal signal and its spectrum in the same figure in two separate graphs (use subplot function).

1.1 Analog low-pass filtering
The goal here is to suppress the three last notes of the signal. The retained specifications are (Rp = 3 dB, Rs = 40 dB, fc = 420 Hz et ∆f = 100 Hz). In order to compare the four types of filters (Butterworth, Chabyshev de type 1 ou 2 et Cauer).

1. Calculate using Matlab the required ordre of each of the four filters. Which one is more advantageous ?
2. Draw on the same graph the magnitude frequency responses (expressed in dB) of the four filter quatre filtres (use the function legend). Give comments regarding the difference between the filter an explain if the required specifications are satisfied.
3. Apply the four filters to the audio signal. Listen the filter signal. Are all the notes suppressed ?
4. For each filter, draw on the same figure the filtered signal and its magnitude spectrum (use the function subplot). The interval of plotting should be limited to 1000 frequencies between 250 Hz and 550 Hz. Comment the results. Which results are more satisfying ?
5. Do the same analysis and comment the results when taking the following specifications :
a) reduction of the transition band width : (Rp = 3 dB and ∆f = 20 Hz).
b) re ́duction of the ripples in the pass-band : (Rp = 1 dB and ∆f = 20 Hz).

# Rejection of one note by digital filtering
The goal here is to eliminate a single note from the audio signal using IIR and FIR filters.
1. Calculate using Matlab the digital filter of type Chebyshev 2 tha allows to eliminate e ́limine the note FA from the signal. We will retain the following specifications (fcl = 340 Hz, fch = 360 Hz, ∆fl = ∆fh = 10Hz,Rp =1dB,Rs =40dB)
2. Draw on the same figure the filtered signal and its magnitude spectrum. We will restrict the frequencies to 1000 equidistant values between 250 Hz and 550 Hz,
3. Do the same operations, as in the previous question, using a RIF filter having the lowest length (Kaiser window).
4. Draw on the same figure the magnitude frequency responses (in decibels) of the RII et FIR filters. Give somme comment about the differences.
