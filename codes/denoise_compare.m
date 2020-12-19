clc; clear all; close all;
tic
%% step 1: load data from station
station = {'lbba','lbse','lbno','lbsw'};
%e stands for east, n stands for north, z stands for vertical
dir = {'e','n','z'};
path = 'C:\Users\think\Desktop\UCSDWC\2020fall\251c\ConvertedData\ConvertedData'; 
date.year = '2020'; date.month = '04'; date.day = '04'; date.hour = '01';
fs=100;start=321000;secs=60;
truncate= [start+1 start+secs*fs];%地震信号
%truncate= [1 40000];
[input,output] = load_data(path, station, dir, date, truncate);
toc
%lbba_e_t=input(:,1);
lbba_e_t=output(:,1);
n=length(input(:,1)); %点数
t=[1:n]/fs;
% figure(1); plot(t,output(:,2));

%% fft 频谱
lbba_e_freq = fft(input(:,1));
fq = (0:n-1)*(fs/n);
lbba_e_power = abs(lbba_e_freq).^2/n;
figure();
plot(fq(1:floor(n/2)),lbba_e_power(1:floor(n/2)))
xlabel('Frequency');ylabel('Power');title('original signal')
%% wavelet spectrogram
figure(3);
cwt(input(:,1),fs);title('cwt for original sig')
%% wavelet denoise 1
lbba_e_den = wdenoise(lbba_e_t,10,'wavelet','db6');
%lbba_e_den = wdenoise(lbba_e_t,'DenoisingMethod','BlockJS');
% cleansym = wdenoise(fNoisy,9,'Wavelet','sym4');
% cleandb = wdenoise(fNoisy,9,'Wavelet','db1');
figure();
subplot(211);
plot(t,lbba_e_t);title('原始信号');
xlabel('time(s)');ylabel('amplitude');
subplot(212);
plot(t,lbba_e_den);title('wave1消噪后的信号');
%h1(2).LineWidth = 2;
xlabel('time(s)');ylabel('amplitude');
figure(4);
cwt(lbba_e_den,fs);title('wave1 denoise cwt')
%% wavelet denoise 2
lev = 10;
lbba_e_den2 = wden(lbba_e_t,'minimaxi','s','mln',lev,'sym4');
figure();
subplot(211);
plot(lbba_e_t);title('原始信号');
xlabel('time(s)');ylabel('amplitude');
subplot(212);
plot(lbba_e_den2);title('wave2消噪后的信号');
xlabel('time(s)');ylabel('amplitude');
figure();
cwt(lbba_e_den2,fs);title('wave2 denoise cwt')

%% STFT denoise
[s,f,t_stft]= stft(lbba_e_t,fs);
figure();
surf(t_stft,f,abs(s),'edgecolor','none'); %有用
% shading flat;
% axis tight;
view(0,90);
% colormap default
% hcol = colorbar;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
ylim([0 fs/2]);
%set(gca, 'FontName', 'Times New Roman')
% set(hcol, 'FontName', 'Times New Roman', 'FontSize', 8);
%threshold = 0.03;  %改
sigma=median((median(abs(s))))/0.6745;
threshold=sigma*sqrt(2*log(7000)); % universal threshold

s(abs(s)<threshold) = 0;
[x_denoise, t_denoise] = istft(s,fs);
figure();
subplot(211);
plot(t, lbba_e_t);
title('Original Signal');
xlabel('Time (sec)');ylabel('Accleration (m/sec^2)');

% hold on;
subplot(212);
plot(t_denoise, x_denoise);
title('Denoised signal(STFT)');
%ylim([min(lbba_e_t)*1.1 max(lbba_e_t)*1.1])
xlabel('Time (sec)');ylabel('Accleration (m/sec^2)');

%set(gca, 'FontName', 'Times New Roman')
