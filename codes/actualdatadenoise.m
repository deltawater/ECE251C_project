clc; clear all; close all;
% This code is to RUN DSI for Geisel Library
%fs=100 %input=lbba(basement), %output:other three
tic
%% step 1: load data from station
station = {'lbba','lbse','lbno','lbsw'};
%e stands for east, n stands for north, z stands for vertical
dir = {'e','n','z'};
path = 'C:\Users\think\Desktop\UCSDWC\2020fall\251c\ConvertedData\ConvertedData'; 
date.year = '2020'; date.month = '04'; date.day = '04'; date.hour = '01';
fs=100;Fs=fs;
start=321000;secs=60;
truncate= [start+1 start+secs*fs];%地震信号
%truncate= [1 40000];
[input,output] = load_data(path, station, dir, date, truncate);
toc
lbba_e=input(:,1);lbba_n=input(:,2);lbba_z=input(:,3);
lbse_e=output(:,1);lbse_n=output(:,2);lbse_z=output(:,3);
lbno_e=output(:,4);lbno_n=output(:,5);lbno_z=output(:,6);
lbsw_e=output(:,7);lbsw_n=output(:,8);lbsw_z=output(:,9);
n=length(input(:,1)); %点数
t=[1:n]/fs;
% figure(1); plot(t,output(:,2));
xblocks_n=lbsw_e;
N=5; %分解层数, snr 最高
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


% blocksoft_snrt=snr(xblocks,xblocks-blocksoft);
% %blocksoft_rmse=sqrt(mse(xblocks-blocksoft));
% blocksoft_rmse=sqrt(sum((xblocks-blocksoft).^2)/1024);
%% denoise with hard
ds_h=hard(D,T);
cs_h=C_block;
for i=1:N
    cs_h(num(i)+1:num(i+1))=ds_h(N-i+1,1:L_block(i+1));
end
blockhard=waverec(cs_h,L_block,'db4');

% blockhard_snrt=snr(xblocks,xblocks-blockhard);
% blockhard_rmse=sqrt(mse(xblocks-blockhard));
% %blockhard_rmse=sqrt(sum((xblocks-blockhard).^2)/1024);

%% denoise with erf
ds_e=erff(D,T,500);
cs_e=C_block;
for i=1:N
    cs_e(num(i)+1:num(i+1))=ds_e(N-i+1,1:L_block(i+1));
end
blockerf=waverec(cs_e,L_block,'db4');

% blockerf_snrt=snr(xblocks,xblocks-blockerf);
% blockerf_rmse=sqrt(mse(xblocks-blockerf));


%% garr
ds_g=garrote(D,T);
cs_g=C_block;
for i=1:N
    cs_g(num(i)+1:num(i+1))=ds_g(N-i+1,1:L_block(i+1));
end
blockgarr=waverec(cs_g,L_block,'db4');


%% show
figure();
title('original signal');
subplot(221)
plot(xblocks_n,'color',[0.7,0.7,0.7]);hold on;plot(blocksoft);
legend('orignal signal','soft-thresh signal');grid on;
title('soft-thresh signal');xlabel('Samples');ylabel('Amplitute');
subplot(222)
plot(xblocks_n,'color',[0.7,0.7,0.7]);hold on;plot(blockhard);
legend('orignal signal','hard-thresh signal');grid on;
title('hard-thresh signal');xlabel('Samples');ylabel('Amplitute');
subplot(223)
plot(xblocks_n,'color',[0.7,0.7,0.7]);hold on;plot(blockgarr);
legend('orignal signal','garrote-thresh signal');grid on;
title('garrote-thresh signal');xlabel('Samples');ylabel('Amplitute');
subplot(224)
plot(xblocks_n,'color',[0.7,0.7,0.7]);hold on;plot(blockerf);
legend('orignal signal','erf-thresh signal');grid on;%beita=500
title('erf-thresh signal');xlabel('Samples');ylabel('Amplitute');

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

