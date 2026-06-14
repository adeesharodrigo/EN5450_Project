%Index number is 259139E
%Initial Parameters

%A=1 B=3 C=9
A=1;
B=3;
C=9;
D=1;
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

%Normalizing frequencies (discrete time frequencies)
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
fprintf('Deviation in passband = %.5f\n', delta_p);
fprintf('Deviation in stopband = %.5f\n', delta_s);
fprintf('Delta = %.4f\n', delta);



% Prewarped frequencies are indicated with pw_x
Nfft=4096;
pw_w_p1=(2/T)*tan(wp1/2);
pw_w_p2=(2/T)*tan(wp2/2);
pw_w_a1=(2/T)*tan(wa1/2);
pw_w_a2=(2/T)*tan(wa2/2);

fprintf('prewarped Lower passband edge = %.4f\n', pw_w_p1);
fprintf('prewarped upper passband edge = %.4f\n', pw_w_p2);
fprintf('prewarped lower stopband edge = %.4f\n', pw_w_a1);
fprintf('prewarped upper stopband edge = %.4f\n', pw_w_a2);


%Finding minimum order analog chebyshev filter
% cheb1ord: calculates the minimum order of an analog or a digital Chebyshev filter required to
%satisfy given specifications
[N_iir, Wn_iir] = cheb1ord([pw_w_p1 pw_w_p2],[pw_w_a1 pw_w_a2], A_p, A_a, 's');    
fprintf('Chebyshev Analog Order N = %d\n', N_iir); 

%Chebyshev bandpass filter
%cheby1: calculates the coefficients of the transfer functions of analog or IIR digital Chebyshev
%filters
%Get Numerator and denominator of H(s) (Analog transfer function)
[Numerator, Denominator] = cheby1(N_iir, A_p, [pw_w_p1 pw_w_p2], 'bandpass', 's');

%Get bilinear tranformation
%bilinear: calculates the coefficients of the transfer function of an IIR digital filter for a given
%transfer function of an analog filter
%Consider numerator of transfer function is B_n,denominator as B_d
[B_n, B_d]=bilinear(Numerator,Denominator,s);

fprintf('Digital IIR filter order = %d\n', length(B_n)-1);  % Order = length of A - 1

%Part 3a: Tabulating
fprintf('\n--- Numerator coefficients B (b0 … bN) ---\n');
fprintf('  b%d = %.8e\n', [(0:length(B_n)-1); B_n(:)']);  % Print each numerator coefficient b_k
fprintf('\n--- Denominator coefficients A (a0 … aN) ---\n');
fprintf('  a%d = %.8e\n', [(0:length(B_d)-1); B_d(:)']);  % Print each denominator coefficient a_k

disp('Coefficient table (IIR – Chebyshev Type I):');
fprintf('%-6s  %-20s  %-20s\n', 'Index', 'Numerator (b)', 'Denominator (a)');  % Table header
max_len = max(length(B_n), length(B_d));                   % Number of rows = longer of B or A
B_pad = [B_n, zeros(1, max_len - length(B_n))];            % Zero-pad B if shorter than A
A_pad = [B_d, zeros(1, max_len - length(B_d))];            % Zero-pad A if shorter than B
for k = 1:max_len
    fprintf('%-6d  %-20.8e  %-20.8e\n', k-1, B_pad(k), A_pad(k));  % Print row k: index, b_k, a_k
end


%  Part 3b: Plot magnitude response from -pi to pi --
[H_iir, w_iir] = freqz(B_n, B_d, Nfft, 'whole');  % Compute digital IIR frequency response over [0, 2pi)
w_iir_shifted  = w_iir - pi;                       % Shift frequency axis to [-pi, pi)
H_iir_shifted  = fftshift(H_iir);                  % Rearrange H(z) to match shifted axis
figure('Name', 'Part 3b – IIR Magnitude Response');                  
plot(w_iir_shifted/pi, 20*log10(abs(H_iir_shifted)));                 % Plot magnitude in dB
xlabel('\omega / \pi  (rad/sample)');                                 % x-axis label
ylabel('|H(e^{j\omega})| (dB)');                                      % y-axis label
title('Part 3b: IIR Chebyshev Magnitude Response  (–\pi to \pi)');  
grid on;                                                              
xlim([-1 1]);                                                         
yline(-A_p, 'r--', 'Passband limit');   % Dashed red line at passband ripple limit
yline(-A_a, 'b--', 'Stopband limit');  % Dashed blue line at stopband attenuation limit

run('C:\Users\adees\Downloads\DSP Project\task1') ;
% Part 3c: Zoom into the passband --
figure('Name', 'Part 3c – IIR Passband');                      % New figure
plot(w_iir(w_pb_idx), 20*log10(abs(H_iir(w_pb_idx))));         % Plot only passband portion of H(z)
xlabel('\omega (rad/sample)');                                  % x-axis label
ylabel('|H(e^{j\omega})| (dB)');                               % y-axis label
title('Part 3c: IIR Chebyshev Passband Detail');                % Title
grid on;                                                        % Enable grid
xline(wp1, 'r--');  % Vertical line at lower passband edge
xline(wp2, 'r--');  % Vertical line at upper passband edge