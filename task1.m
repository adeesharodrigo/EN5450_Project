%Index number is 259139E
%Initial Parameters

%A=1 B=3 C=9
A=1;
B=3;
C=9;
%Maximum Passband Ripple
A_p=0.03+(0.01*A); %in dB
%Minimum Stopband attenuation
A_a=40+B; %in dB
%Lower passband edge (W_p1)
W_p1=(C*100)+400; %rad/s
%Upper passband edge
W_p2=(C*100)+900; %rad/s
%Lower stopband edge (W_a1)
W_a1=(C*100)+200; %rad/s
%Upperstopband edge (W_a2)
W_a2=(C*100)+1100; %rad/s
 %Sampling Frequency in rad/s
W_s= 2*((C*100)+1600) ;
% Finding Nyquist Frequency
Nyquist=W_s/2; %in rad/s

s=Nyquist/pi; % Sampling frequency in Hz
T=1/s; %Sampling period

fprintf(' Lower passband edge = %.4f\n', W_p1);
fprintf(' upper passband edge = %.4f\n', W_p2);
fprintf(' lower stopband edge = %.4f\n', W_a1);
fprintf(' upper stopband edge = %.4f\n', W_a2);
fprintf(' Nyquist frequency = %.4f\n', Nyquist);
fprintf(' Sampling frequency in Hz = %.4f\n', s);

%Normalizing frequencies
wp1=W_p1*T;
wp2=W_p2*T;
wa1=W_a1*T;
wa2=W_a2*T;

fprintf('Maximum Passband Ripple = %.2f\n', A_p);
fprintf('Minimum Stopband attenuation = %.2f\n', A_a);
fprintf('Normalized Lower passband edge = %.4f\n', wp1);
fprintf('Normalized upper passband edge = %.4f\n', wp2);
fprintf('Normalized lower stopband edge = %.4f\n', wa1);
fprintf('Normalized upper stopband edge = %.4f\n', wa2);


%dB to linear conversion
delta_p=(10^(A_p/20)-1)/(10^(A_p/20)+1);  %Maximum passband deviation (linear)
delta_s=(10^(-A_a/20));                    %Maximum stopband deviation (lieanr)
delta=min(delta_s,delta_p);%worst case deviation for kaiser rod
fprintf('Deviation in passband (delta_p) = %.2f\n', delta_p);
fprintf('Deviation in stopband (delta_s)= %.2f\n', delta_s);
fprintf('Delta = %.2f\n', delta);

%% Kaizer filter

f_edges = [wa1 wp1 wp2 wa2]/pi; % defining band edges
mags = [0 1 0]; %desired magnitudes in stop band, passband and stop band
devs = [delta_s delta_p delta_s]; %Maximum allowed deviations in each band

%returns a filter order n, normalized frequency band edges Wn, and a shape factor beta that specify a Kaiser window for use with the fir1 function
[N_kaiser, Wn_kaiser, beta_kaiser, ftype_kaiser] = kaiserord(f_edges, mags, devs);
%Validate whether N_kaiser is even
if mod(N_kaiser, 2) ~= 0
    N_kaiser = N_kaiser + 1;
end
fprintf('Order of Kaizer filter = %.2f\n', N_kaiser);
fprintf('Shape factor = %.2f\n',beta_kaiser);

%Filter
hh = fir1(N_kaiser,Wn_kaiser,ftype_kaiser,kaiser(N_kaiser+1,beta_kaiser),"noscale");

% [H,f] = freqz(hh,1,1024,s);
% plot(f,abs(H))
% grid


% a) Plotig impulse response
figure('Name','Part 1a – Kaiser Impulse Response');
n_vec = 0:N_kaiser;  %Sample Index vector
stem(n_vec, hh, 'filled', 'MarkerSize', 3); %ploting impulse response
xlabel('n (samples)'); ylabel('h[n]');
title(sprintf('Part 1a: Kaiser FIR Impulse Response  (N=%d)', N_kaiser));
grid on;

% b) Magnitude response

Nfft = 4096;  %number of FFT points
[H_k, w_k] = freqz(hh, 1, Nfft, 'whole');  %Frequencies form 0 to 2 pi
w_k_shifted = w_k - pi;              % shift 0 to 2 pi to -pi to pi
H_k_shifted = fftshift(H_k);
[H_k, w_k] = freqz(hh, 1, Nfft, 'whole');


figure('Name','Part 1b – Kaiser Magnitude Response');
plot(w_k_shifted/pi, 20*log10(abs(H_k_shifted)));
xlabel('\omega / \pi  (rad/sample)'); ylabel('|H(e^{j\omega})| (dB)');
title('Part 1b: Kaiser FIR Magnitude Response  (–\pi to \pi)');
grid on; xlim([-1 1]);
yline(-A_p,  'r--', 'Passband limit', 'LabelHorizontalAlignment','left');
yline(-A_a, 'b--', 'Stopband limit',  'LabelHorizontalAlignment','left');


% c) Passband
figure('Name','Part 1c – Kaiser Passband');
w_pb_idx = (w_k >= wp1) & (w_k <= wp2);% Boolean index for passband frequencies
plot(w_k(w_pb_idx), 20*log10(abs(H_k(w_pb_idx))));
xlabel('\omega (rad/sample)'); ylabel('|H(e^{j\omega})| (dB)');
title('Part 1c: Kaiser FIR Passband Detail');
grid on;
xline(wp1,'r--'); xline(wp2,'r--');

