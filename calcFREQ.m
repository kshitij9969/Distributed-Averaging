function freq = calcFREQ(s,h,t)
%%This function computes the fast fourier transform of the switching
%%function s.
d=s;
Ts = h;                                                     % Sampling Interval (sec)
Fs = 1/Ts;                                                  % Sampling Frequency (Hz)
Fn = Fs/2;                                                  % Nyquist Frequency (Hz)
L = length(t);                                             % Vector Length
FTs = fft(d-mean(d))/L;                                    % Fourier Transform (Subtract d-c Offset)
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
Iv = 1:length(Fv);                                          % Index Vector
[pks1,frqs1] = findpeaks(abs(FTs(Iv))*2, Fv, 'MinPeakHeight',0.4);


% fund_freq=linspace(0,frqs1*1.5);
freq=frqs1((pks1==max(pks1)))
figure(1)
plot(Fv, (abs(FTs(Iv))*2))
xlim([0 freq*1.5]);
hold on
plot(frqs1, pks1, '^b')
grid
axis([0  freq*1.5    ylim])


end