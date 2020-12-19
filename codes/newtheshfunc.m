% new thresholding function erf/tanh
clc; clear all; close all;
% use 4 standard signal to evaluate (from SNR and RMSE)
% (a) Blocks signal; (b) Bumps signal; (c) Heavy sine signal; (d) Doppler signal.
% and simulated vibration signal and true data(geisel library)
%three traditional thresholding function:‘soft’, ‘hard’, ‘garrote’
%% thresholding functions
% x=-5:1/100:5;T=1;beita=1;
% figure();
% % ytanh=tanhf(x,T,beita);
% % plot(x,ytanh);hold on;
% yerf=erff(x,T,beita);
% plot(x,yerf,'linewidth',2);hold on;
% ygarr=garrote(x,T);
% plot(x,ygarr,'linewidth',2);hold on;
% ysoft=soft(x,T);yhard=hard(x,T);
% %figure();
% plot(x,ysoft,'linewidth',2);grid on;
% hold on;
% plot(x,yhard,'linewidth',2);
% % legend('tanh','erf','garrote','soft','hard');
% legend('erf','garrote','soft','hard');
% title('Comparison of different thresholding functions');
%% generate 4 test signals (a) Blocks signal; (b) Bumps signal; (c) Heavy sine signal; (d) Doppler signal.
xblocks=wnoise('blocks',10); xbumps=wnoise('bumps',10);
xhsine=wnoise('heavy sine',10); xdoppler=10*wnoise('doppler',10);
% figure();subplot(221);plot(xblocks);title('blocks signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(222);plot(xbumps);title('bumps signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(223);plot(xhsine);title('heavy sine signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(224);plot(xdoppler);title('doppler signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;

%% add white noise to 4 standard signals
SNR=2;
xblocks_n = awgn(xblocks,SNR,'measured');xbumps_n = awgn(xbumps,SNR,'measured');
xhsine_n = awgn(xhsine,SNR,'measured');xdoppler_n = awgn(xdoppler,SNR,'measured');
% figure();subplot(221);plot(xblocks_n);title('noisy blocks signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(222);plot(xbumps_n);title('noisy bumps signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(223);plot(xhsine_n);title('noisy heavy sine signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
% subplot(224);plot(xdoppler_n);title('noisy doppler signal');xlim([0 1024]);xlabel('Samples'); ylabel('Amplitude');grid on;
xblocks=xdoppler;xblocks_n=xdoppler_n;
bori_rmse=sqrt(mse(xblocks-xblocks_n));
%% denoise with various thresh functions
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
% T2=thselect(xblocks_n,'minimaxi');
% T=thselect(D(1,:),'minimaxi');
%% soft
ds_s=soft(D,T);
cs_s=C_block;
for i=1:N
    cs_s(num(i)+1:num(i+1))=ds_s(N-i+1,1:L_block(i+1));
end
blocksoft=waverec(cs_s,L_block,'db4');

%blocksoft_snr=msnr(xblocks,xblocks_n);
blocksoft_snrt(N)=snr(xblocks,xblocks-blocksoft);
blocksoft_rmse(N)=sqrt(mse((xblocks-blocksoft)));
% blocksoft_rmset=sqrt(sum((xblocks-blocksoft).^2)/1024);
%% denoise with hard
ds_h=hard(D,T);
cs_h=C_block;
for i=1:N
    cs_h(num(i)+1:num(i+1))=ds_h(N-i+1,1:L_block(i+1));
end
blockhard=waverec(cs_h,L_block,'db4');

blockhard_snrt(N)=snr(xblocks,xblocks-blockhard);
blockhard_rmse(N)=sqrt(mse((xblocks-blockhard)));


%% denoise with erf
ds_e=erff(D,T,1);
cs_e=C_block;
for i=1:N
    cs_e(num(i)+1:num(i+1))=ds_e(N-i+1,1:L_block(i+1));
end
blockerf=waverec(cs_e,L_block,'db4');

blockerf_snrt(N)=snr(xblocks,xblocks-blockerf);
blockerf_rmse(N)=sqrt(mse((xblocks-blockerf)));

%% denoise with garrote
ds_g=garrote(D,T);
cs_g=C_block;
for i=1:N
    cs_g(num(i)+1:num(i+1))=ds_g(N-i+1,1:L_block(i+1));
end
blockgarr=waverec(cs_g,L_block,'db4');
blockgarr_snrt(N)=snr(xblocks,xblocks-blockgarr);
blockgarr_rmse(N)=sqrt(mse((xblocks-blockgarr)));

% end
%% show
% 加上noised和原始的
figure();
plot(blocksoft,'linewidth',1);
plot(xblocks_n,'color',[0.75,0.75,0.75]);hold on;
plot(xblocks,'linewidth',1.2);hold on;
plot(blocksoft,'linewidth',1.2);hold on;
%plot other3 thresh functions
plot(blockhard,'linewidth',1.2);hold on;
plot(blockgarr,'linewidth',1.2);hold on;
plot(blockerf,'linewidth',1.2);
legend('noisy signal','original signal','soft-thresh signal','hard-thresh signal','garrote-thresh signal','erf-thresh signal');
xlabel('Samples');ylabel('Amplitute');xlim([0 1024]);
title('Comparison of different threshold function denoise result');
% subplot(222)
% plot(blockhard);
% title('blocks hard-thresh');
% subplot(223)
% plot(blockerf);
% title('blocks erf-thresh');
% subplot(224)
% plot(blockgarr);
% title('blocks garrote-thresh');

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

