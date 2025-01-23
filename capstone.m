clear all;
clc;
rx = adi.Pluto.Rx('uri','ip:Pluto');
tx = adi.Pluto.Tx('uri','ip:Pluto');
%%
samp_rate = 2e6;   % # must be <=30.72 MHz if both channels are enabled
NumSamples = 2^12;
rx_lo = 2.3e9;
rx_mode = "manual"; % # can be "manual" or "slow_attack"
rx_gain0 = 40;
rx_gain1 = 40;
tx_lo = rx_lo;
tx_gain = 0;
fc0 = (200e3);
num_scans = 5;
phase_cal = 160;
Plot_Compass = true;



d_wavelength = 0.5           ;       % distance between elements as a fraction of wavelength.  This is normally 0.5
wavelength = 3e8/rx_lo      ;        % wavelength of the RF carrier
d = d_wavelength*wavelength;    
signal_start = int16(NumSamples*(samp_rate/2+fc0/2)/samp_rate);
signal_end = int16(NumSamples*(samp_rate/2+fc0*2)/samp_rate)  ;                      % distance between elements in meters

%%

rx.EnabledChannels = [1,2]; %0 1 olabilir
rx.SamplingRate = samp_rate;
rx.RFBandwidth = fc0*3;
rx.CenterFrequency= rx_lo;
rx.GainControlModeChannel0 = "manual";
rx.GainControlModeChannel1 = "manual";
rx.GainChannel0 = rx_gain0;
rx.GainChannel1 = rx_gain1;
rx.SamplesPerFrame =NumSamples;
rx.kernelBuffersCount = 1;



tx.RFBandwidth = fc0*3;
tx.CenterFrequency = tx_lo;
tx.EnableCyclicBuffers = true;
tx.AttenuationChannel0 = tx_gain;

% tx.SamplesPerFrame = NumSamples;
tx.SamplingRate = samp_rate;
fs =(rx.SamplingRate);
N = 2^16;
ts = 1 / (fs);
t = 0:ts:(N * ts - ts)*2;

i0 = cos(2 * pi * t * fc0) * 2^14;
q0 = sin(2 * pi * t * fc0) * 2^14;
iq0 = i0 + 1j * q0;
iq0 = iq0.';
amplitude = 2^15; 
frequency = 0.17e6;
swv1 = dsp.SineWave(amplitude, frequency);
swv1.ComplexOutput = true;
swv1.SamplesPerFrame = 1e4*10;
swv1.SampleRate = 3e6;
y = swv1();
tx.step([iq0;iq0]);
%# Send Tx data.


%%
%tx.release();
  %
data = rx();
figure
subplot(3,1,1);
plot((abs(double(data(:,1))))),
hold on


plot((abs(double(data(:,2)))))
fs = rx.SamplingRate;


L = length(data(:,1));
f=linspace(-fs/2,fs/2,L);
Y1=fftshift(fft(data(:,1),L)/L);
subplot(3,1,2);
plot(f, abs(Y1));
max(abs(Y1))

title("Channel 1 (Default)");
Y2=fftshift(fft(data(:,2),L)/L);

subplot(3,1,3);
plot(f, abs(Y2));

title("Channel 2 (uFL)")

xlim([-1e6 1e6]);
xlabel("Frequency (Hz)")

% max(abs(Y2))
% rx.release()
% tx.release()
%tx.release();


%%
for i=1:10
    temp_data_for_calibration = rx();
end
steer_angles = zeros([1 num_scans]);
peak_delays = zeros([1 num_scans]);
for i=1:num_scans
    data = rx();
    rx_0 = data(:,1);
    rx_1 = data(:,2);
    peak_sum = zeros(1 , 181); %asagıdaki phase delay loopunun incrementine baglı sizeı var.
    temporary_index_for_phase_loop = 1;
    delay_phases = -180:2:180;
    for phase_delay=-180:2:180
        delay_rx1 = rx_1.*exp(1i*deg2rad(phase_delay + phase_cal));
        delayed_sum = dbfs(rx_0+delay_rx1);

        peak_sum(1,temporary_index_for_phase_loop) = max(delayed_sum(signal_start:signal_end,1));
        temporary_index_for_phase_loop  = temporary_index_for_phase_loop +1;
    end
    peak_dbfs = max(peak_sum);
    peak_delay_index = find(peak_sum==peak_dbfs);
    peak_delay = delay_phases(peak_delay_index(1));
    steer_angle = int16(calcTheta(peak_delay,d,rx_lo));
    figure;

    plot(delay_phases, peak_sum);
    
    % line([peak_delay, peak_delay], ylim, 'Color', 'g', 'LineStyle', ':');
    xline(peak_delay,'Color', 'r', 'LineStyle', ':')
    % text(-180, -26, sprintf('Peak signal occurs with phase shift = %.1f deg', round(peak_delay, 1)));
    % text(-180, -28, sprintf('If d=%d mm, then steering angle = %.1f deg', int32(d*1000), steer_angle));
    ylim([-80 0]);
    xlim([-180 180]);
    xlabel('phase shift [deg]');
    ylabel('Rx0 + Rx1 [dBfs]');
    peak_delays(i) = peak_delay;
    steer_angles(i) = steer_angle;
end

steer_angles
peak_delays

    