
% load('BER_MMSE.mat');
% load('BER_PE_3.mat');
% load('BER_LPE_3.mat');
% %%
% % 图一
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSE detect','L=3 PE detect','L=3 LPE detect');
% %%
% 图二
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% hold on;
% load('BER_PE_2.mat');
% load('BER_LPE_2.mat');
% semilogy(SNR,BER_PE,'Color','red','LineStyle','--','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','--','Marker','*');
% load('BER_PE_1.mat');
% load('BER_LPE_1.mat');
% semilogy(SNR,BER_PE,'Color','red','LineStyle','-.','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-.','Marker','*');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSE detect','L=3 PE detect','L=3 LPE detect','L=2 PE detect','L=2 LPE detect','L=1 PE detect','L=1 LPE detect');

%%
% %图三
% semilogy(SNR,BER_LPE,'Color','red','LineStyle','-','Marker','*');
% hold on;
% load('BER_LPE_2.mat');
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','o');
% hold on;
% load('BER_LPE_1.mat');
% semilogy(SNR,BER_LPE,'Color','blue','LineStyle','-','Marker','+');
% legend('L=3 LPE detect','L=2 LPE detect','L=1 LPE detect');


%%
% semilogy(SNR,BER_BDMMSE,'Color','green','LineStyle','-','Marker','^');
% hold on;
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSEdetect');

%%
figure(1);
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR,BER_PE(:,1),'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,2),'Color','red','LineStyle','--','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,3),'Color','red','LineStyle',':','Marker','+');
xlabel('SNR(in dB)');
ylabel('BER');
legend('MMSEdetect','L=1 PEdetect','L=2 PEdetect','L=3 PEdetect');

%%
figure(2);
semilogy(SNR(1:6),BER_MMSE(1:6),'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR(1:6),BER_LPE(1:6,1),'Color','black','LineStyle','-','Marker','s');
hold on;
semilogy(SNR(1:6),BER_LPE(1:6,2),'Color','black','LineStyle','--','Marker','s');
hold on;
semilogy(SNR(1:6),BER_LPE(1:6,3),'Color','black','LineStyle',':','Marker','s');
hold on;
semilogy(SNR(1:6),BER_PE(1:6,1),'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR(1:6),BER_PE(1:6,2),'Color','red','LineStyle','--','Marker','+');
hold on;
semilogy(SNR(1:6),BER_PE(1:6,3),'Color','red','LineStyle',':','Marker','+');
xlabel('SNR(in dB)');
ylabel('BER');
legend('MMSEdetect','L=1 LPEdetect','L=2 LPEdetect','L=3 LPEdetect','L=1 PEdetect','L=2 PEdetect','L=3 PEdetect');
%%
figure(3);
semilogy(SNR,BER_LPE(:,1),'Color','black','LineStyle','-','Marker','s');
hold on;
semilogy(SNR,BER_LPE(:,2),'Color','black','LineStyle','--','Marker','s');
hold on;
semilogy(SNR,BER_LPE(:,3),'Color','black','LineStyle',':','Marker','s');
hold on;
semilogy(SNR,BER_PE(:,1),'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,2),'Color','red','LineStyle','--','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,3),'Color','red','LineStyle',':','Marker','+');
xlabel('SNR(in dB)');
ylabel('BER');
legend('L=1 LPEdetect','L=2 LPEdetect','L=3 LPEdetect','L=1 PEdetect','L=2 PEdetect','L=3 PEdetect');

%%
% load('BER_MMSE1.mat')
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% load('BER_MMSE1.mat')
% semilogy(SNR,BER_MMSE,'Color','red','LineStyle','-','Marker','^');
% hold on;
% xlabel('SNR(in dB)');
% ylabel('BER');
% legend('MMSEdetect','QR分解 MMSEdetect');