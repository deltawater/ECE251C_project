clc; clear all; close all;
% This code is to RUN DSI for Geisel Library
%fs=100 %input=lbba(basement), %output:other three
%% step 1: load data from station
station = {'lbba','lbse','lbno','lbsw'};
%e stands for east, n stands for north, z stands for vertical
dir = {'e','n','z'};
path = 'C:\Users\think\Desktop\UCSDWC\2020fall\251c\simulited_data.xlsx'; 
fs=100;Fs=fs;
start=321000;secs=60;
truncate= [start+1 start+secs*fs];%地震信号

fileID = fopen(path,'r');
formatSpec = '%f';
sdata=xlsread(path);
lbse_e=sdata(4:6004,2)/384.09;
% lbse_e=sdata(14005:14005+6000,3)/384.09;
% lbse_e=sdata(28006:28006+6000,4)/384.09;
% lbse_e=output(:,1);
n=length(lbse_e); %点数
t=[1:n]/fs;
% figure(1); plot(t,output(:,2));
xblocks=lbse_e;
xblocks_n=awgn(lbse_e,12,'measured');
for N=1:4 %分解层数, snr 最高
[C_block,L_block]=wavedec(xblocks_n,N,'db4');
num=1:N+1;
num(1)=L_block(1);
for i=2:N+1
    num(i)=num(i-1)+L_block(i);
end
D=zeros(N,L_block(N+1));
for i=1:N
    D(i,1:L_block(N-i+2))=detcoef(C_block,L_block,i);
end
sigma=(median(abs(D(1,:))))/0.6745;
T=sigma*sqrt(2*log(L_block(N+1))); % universal threshold
%T=thselect(D(1,:),'heursure');
%% soft
ds_s=soft(D,T);
cs_s=C_block;
for i=1:N
    cs_s(num(i)+1:num(i+1))=ds_s(N-i+1,1:L_block(i+1));
end
blocksoft=waverec(cs_s,L_block,'db4');


blocksoft_snrt(N)=snr(xblocks,xblocks-blocksoft);
blocksoft_rmse(N)=sqrt(mse(xblocks-blocksoft));
% blocksoft_rmse=sqrt(sum((xblocks-blocksoft).^2)/1024);
%% denoise with hard
ds_h=hard(D,T);
cs_h=C_block;
for i=1:N
    cs_h(num(i)+1:num(i+1))=ds_h(N-i+1,1:L_block(i+1));
end
blockhard=waverec(cs_h,L_block,'db4');

blockhard_snrt(N)=snr(xblocks,xblocks-blockhard);
blockhard_rmse(N)=sqrt(mse(xblocks-blockhard));
% %blockhard_rmse=sqrt(sum((xblocks-blockhard).^2)/1024);

%% denoise with erf
ds_e=erff(D,T,200);
cs_e=C_block;
for i=1:N
    cs_e(num(i)+1:num(i+1))=ds_e(N-i+1,1:L_block(i+1));
end
blockerf=waverec(cs_e,L_block,'db4');

blockerf_snrt(N)=snr(xblocks,xblocks-blockerf);
blockerf_rmse(N)=sqrt(mse(xblocks-blockerf));
%% garrote
ds_g=garrote(D,T);
cs_g=C_block;
for i=1:N
    cs_g(num(i)+1:num(i+1))=ds_g(N-i+1,1:L_block(i+1));
end
blockgarr=waverec(cs_g,L_block,'db4');
blockgarr_snrt(N)=snr(xblocks,xblocks-blockgarr);
blockgarr_rmse(N)=sqrt(mse((xblocks-blockgarr)));


end

%%  STFT denoise

[s,f,t_stft]= stft(xblocks,fs);
% figure();
% surf(t_stft,f,abs(s),'edgecolor','none'); %有用
% view(0,90);
% xlabel('Time (sec)');
% ylabel('Frequency (Hz)');
% ylim([0 fs/2]);
%set(gca, 'FontName', 'Times New Roman')
% set(hcol, 'FontName', 'Times New Roman', 'FontSize', 8);
%threshold = 0.03;  %改
sigma=median((median(abs(s))))/0.6745;
threshold=sigma*sqrt(2*log(6000)); % universal threshold

s(abs(s)<threshold) = 0;
[x_denoise, t_denoise] = istft(s,fs);
figure();
% subplot(211);
% plot(t, xblocks);
% title('Original Signal');
% xlabel('Time (sec)');ylabel('Accleration (m/sec^2)');
% 
% % hold on;
subplot(211);
plot(t, blockhard);
title('denoised signal(wavelet with hard threshold)');
xlabel('Time (sec)');ylabel('Accleration (m/sec^2)');xlim([0 60]);ylim([-0.5 0.5]);
subplot(212);
plot(t_denoise, x_denoise);
title('denoised signal(STFT with hard threshold)');
%ylim([min(lbba_e_t)*1.1 max(lbba_e_t)*1.1])
xlabel('Time (sec)');ylabel('Accleration (m/sec^2)');xlim([0 60]);ylim([-0.5 0.5]);
delen=length(x_denoise);
blockstft_snrt(N)=snr(xblocks(1:delen),xblocks(1:delen)-x_denoise);
blockstft_rmse(N)=sqrt(mse((xblocks(1:delen)-x_denoise)));


%% show
% figure();
% subplot(321)
% plot(xblocks);
% title('original signal');
% subplot(322)
% plot(xblocks_n);
% title('signal with noise');
% subplot(323)
% plot(blocksoft);
% title('denoised signal with soft-thresh');
% subplot(324)
% plot(blockhard);
% title('denoised signal with hard-thresh');
% subplot(326)
% plot(blockerf);
% title('denoised signal with erf-thresh');
% subplot(325)
% plot(blockgarr);
% title('denoised signal with garrote-thresh');

%% functions

function x = soft(b,T)
    x = sign(b).*max(abs(b) - T,0);
end
function x = hard(b,T)
    sel = (abs(b)>T);
    x = b.*sel;
end
function x = garrote(b,T)  
    duan=sign(abs(b)-T);
    duan(duan>0)=1;duan(duan<=0)=0;
    x=duan.*(b-(T^2)./b);
end
function x = tanhf(b,T,beita)
    if (nargin<3)
        beita = 1;
    end
    duan=sign(abs(b)-T);
    duan(duan>0)=1;duan(duan<=0)=0;
    sgn=sign(b);
    %x=sgn.*duan.*(abs(b)-T-T*tanh((abs(b)-T)));
    x=duan.*b.*tanh(beita*(abs(b)-T));
end
function x = erff(b,T,beita) %more parameters??
    if (nargin<3)
        beita = 1;
    end
    duan=sign(abs(b)-T);
    duan(duan>0)=1;duan(duan<=0)=0;
    sgn=sign(b);
    %x=sgn.*duan.*(abs(b)-T-T*tanh((abs(b)-T)));
    x=duan.*b.*erf(beita*(abs(b)-T));
end

