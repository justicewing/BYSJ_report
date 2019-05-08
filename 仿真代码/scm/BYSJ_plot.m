load('BER_PE_3.mat');
load('BER_MMSE.mat');
for k = 1:7
    SNR(k) = -20+k*5;
end

semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
xlabel('SNR');
ylabel('BER');
legend('PEdetect','MMSEdetect');