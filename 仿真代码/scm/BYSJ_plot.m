
% load('BER_MMSE.mat');
% load('BER_PE_3.mat');
% load('BER_LPE_3.mat');
% %%
% % Í¼Ò»
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSE detect','L=3 PE detect','L=3 LPE detect');
% %%
% Í¼¶þ
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
% %Í¼Èý
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
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR,BER_PE(:,1),'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,2),'Color','red','LineStyle','--','Marker','+');
hold on;
semilogy(SNR,BER_PE(:,3),'Color','red','LineStyle',':','Marker','+');
xlabel('SNR');
ylabel('BER');
legend('MMSEdetect','L=1 LPEdetect','L=2 LPEdetect','L=3 LPEdetect');

%%
% semilogy(SNR,BER_PE(:,1),'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_PE(:,2),'Color','red','LineStyle','--','Marker','+');
% hold on;
% semilogy(SNR,BER_PE(:,3),'Color','red','LineStyle',':','Marker','+');
% hold on;
% semilogy(SNR,BER_PE(:,4),'Color','red','LineStyle','-.','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,1),'Color','black','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,2),'Color','black','LineStyle','--','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,3),'Color','black','LineStyle',':','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,4),'Color','black','LineStyle','-.','Marker','+');
% xlabel('SNR');
% ylabel('BER');
% legend('L=1 PEdetect','L=2 PEdetect','L=3 PEdetect','L=4 PEdetect','L=1 LPEdetect','L=2 LPEdetect','L=3 LPEdetect','L=4 LPEdetect');