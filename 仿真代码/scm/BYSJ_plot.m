load('BER_PE_3.mat');
load('BER_MMSE.mat');
load('BER_LPE_3.mat');
for k = 1:7
    SNR(k) = -20+k*5;
end

semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
xlabel('SNR');
ylabel('BER');
legend('PE detect','MMSE detect','LPE detect');