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
fprintf('Deviation in passband = %.5f\n', delta_p);
fprintf('Deviation in stopband = %.5f\n', delta_s);
fprintf('Delta = %.4f\n', delta);


%% Parks-McClellan method

f_edges = [wa1 wp1 wp2 wa2]/pi; % defining band edges
mags = [0 1 0]; %desired magnitudes in stop band, passband and stop band
devs = [delta_s delta_p delta_s]; %Maximum allowed deviations in each band


%[n,fo,ao,w] = firpmord(f,a,dev) returns the order estimate n, normalized frequency band edges fo, frequency band amplitudes ao, and weights w that meet input specifications f, a, and dev.
[N_pm, fo_pm, ao_pm, w_pw]=firpmord(f_edges, mags, devs);
fprintf('Oder of the filter = %.2f\n', N_pm);
%Validate order of filter is even

if mod(N_pm, 2) ~= 0
    N_pm = N_pm + 1;  % Ensure even filter order for linear-phase bandpass FIR
end
% fprintf('edge = %.2f\n', f_edges);


%print calculated order
fprintf('Oder of the filter = %.2f\n', N_pm);

%Compute coefficients of filter (Impulse response - As per first filter
%design lecture)
h_filter=firpm(N_pm, fo_pm, ao_pm, w_pw);

% 2 a)Plot the impulse response ---
figure('Name', 'Parks-McClellan method Impulse Response');
stem(0:N_pm, h_filter, 'filled', 'MarkerSize', 3);        % Stem plot of PM filter impulse response
xlabel('n (samples)');                                    % x-axis label
ylabel('h[n]');                                           % y-axis label
title(sprintf(' Parks-McClellan FIR Impulse Response  (N=%d)', N_pm));  % Title with order
grid on;                          % Enable grid


% 2 b) Plot magnitude response from -pi to pi 
Nfft = 4096; 
[H_pm, w_pm_vec] = freqz(h_filter, 1, Nfft, 'whole');  % Compute frequency response over [0, 2pi)
w_pm_shifted = w_pm_vec - pi;    % Shift frequency axis to [-pi, pi)
H_pm_shifted = fftshift(H_pm); % Rearrange H to match shifted frequency axis
 
figure('Name', 'Parks-McClellan method Magnitude Response');   
plot(w_pm_shifted/pi, 20*log10(abs(H_pm_shifted))); % Magnitude in dB
xlabel('\omega / \pi  (rad/sample)');                                   
ylabel('|H(e^{j\omega})| (dB)');                                        
title(' Parks-McClellan FIR Magnitude Response  (–\pi to \pi)'); 
grid on;                                                                
xlim([-1 1]);                                                            % Show full normalised range
yline(-A_p, 'r--', 'Passband limit');   % Dashed red line at passband ripple limit
yline(-A_a, 'b--', 'Stopband limit');  % Dashed blue line at stopband attenuation limit


% 2 c) Plot  magnitude response for the frequencies in the passband.
figure('Name', 'Parks-McClellan method Passband'); 
w_pb_idx = (w_pm_vec >= wp1) & (w_pm_vec <= wp2);  % Same as the task 1
plot(w_pm_vec(w_pb_idx), 20*log10(abs(H_pm(w_pb_idx))));       % Plot passband only (reuse index from Part 1c)
xlabel('\omega (rad/sample)');                                  
ylabel('|H(e^{j\omega})| (dB)');                              
title('Parks-McClellan FIR Passband ');        
grid on;                                                       
xline(wp1, 'r--');  %  lower passband edge
xline(wp2, 'r--');  %  upper passband edge


