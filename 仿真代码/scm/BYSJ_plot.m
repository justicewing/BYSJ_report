
for k = 1:7
    SNR(k) = -20+k*5;
end

load('BER_MMSE.mat');
load('BER_PE_3.mat');
load('BER_LPE_3.mat');
%%
%ͼһ
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% hold on;

%%
%ͼ��
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
%ͼ��
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% hold on;