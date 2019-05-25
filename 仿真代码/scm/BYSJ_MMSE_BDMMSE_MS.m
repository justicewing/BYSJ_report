run('BYSJ_channelmodel.m');
% 信号生成
sample = rand(NoSamples,Nu,Nt);
sample = round(sample);
signal_t = sample;
signal_t(signal_t(:,:,:)==0) = -1;
% signal_t = ones(NoSamples,Nu);
% for i = 1:NoSamples
%     for j = 1:Nu
%         if(sample(i,j) == 0)
%             signal_t(i,j) = -1;
%         end
%     end
% end
%%
% 接收信号
signal_c = zeros(Nr,NoSamples);
for n_link=1:Nu
    for n_path=1:NoPath
        for n_sample=1:NoSamples
            H_i = H(:,:,n_path,n_sample,n_link);
            H_i = H_i(:,:);
            signal_c(:,n_sample) = signal_c(:,n_sample) + H_i'* squeeze(signal_t(n_sample,n_link,:));
        end
    end
end
    
SNR =zeros(7,1);
BER_MMSE =zeros(7,1);
BER_BDMMSE =zeros(7,1);
for k = 1:7
    SNR(k) = -20+k*2;
    sigma2 = 10^(SNR(k)/10);
    %加噪
    signal_r =awgn(signal_c,SNR(k));
    %%
    % MMSE检测

    signal_MMSE = zeros(NoSamples,Nu*Nt);
    H_k =sum(H,3);
    for n_sample=1:NoSamples
            H_A = H_k(:,:,1,n_sample,1);
            for n_link=2:Nu
                H_A =[H_A;H_k(:,:,1,n_sample,n_link)];
            end
         W = (H_A'*H_A + eye(Nr)/sigma2)\H_A';
         signal_MMSE(n_sample,:) = W' *signal_r(:,n_sample);
    end
    
    
    %多流判决
     rece_MMSE = zeros(NoSamples,Nu,Nt);
     for n_link = 1:Nu
         for n_sample = 1:NoSamples
              for n_trans = 1:Nt
                  rece_MMSE(n_sample,n_link,n_trans) = signal_MMSE(n_sample,Nt*(n_link-1)+n_trans);
              end
          end
      end
    result_MMSE = ones(NoSamples,Nu,Nt);
    for i = 1:NoSamples
        for j = 1:Nu
            for p = 1:Nt
                d_1 = abs(rece_MMSE(i,j,p)-1);
                d_2 = abs(rece_MMSE(i,j,p)+1);
                if(d_1 > d_2)
                    result_MMSE(i,j,p) = 0;
                end
            end
        end
    end
    errortimes = sum(abs(result_MMSE-sample),'all');
    BER_MMSE(k) = sum(errortimes)/(Nu*Nt*NoSamples);
    % BDMMSE检测

    signal_BDMMSE = zeros(NoSamples,Nu*Nt);
    H_k=sum(H,3);
    for n_sample=1:NoSamples
        H_A = H_k(:,:,1,n_sample,1);
        for n_link=2:Nu
            H_A =[H_A;H_k(:,:,1,n_sample,n_link)];
        end
        H_A =H_A*Ut;
        H_a =abs(H_A);
        H_A(H_a(:,:)<0.5) =0;
         W = (H_A'*H_A + eye(Nr)/sigma2)\H_A';
         signal_BDMMSE(n_sample,:) = W' *Ut'*signal_r(:,n_sample);
         
    end
    
    %多流判决
     rece_BDMMSE = zeros(NoSamples,Nu,Nt);
     for n_link = 1:Nu
         for n_sample = 1:NoSamples
              for n_trans = 1:Nt
                  rece_BDMMSE(n_sample,n_link,n_trans) = signal_BDMMSE(n_sample,Nt*(n_link-1)+n_trans);
              end
          end
      end
    result_BDMMSE = ones(NoSamples,Nu,Nt);
    for i = 1:NoSamples
        for j = 1:Nu
            for p = 1:Nt
                d_1 = abs(rece_BDMMSE(i,j,p)-1);
                d_2 = abs(rece_BDMMSE(i,j,p)+1);
                if(d_1 > d_2)
                    result_BDMMSE(i,j,p) = 0;
                end
            end
        end
    end
    errortimes = sum(abs(result_BDMMSE-sample),'all');
    BER_BDMMSE(k) = sum(errortimes)/(Nu*Nt*NoSamples);
    
    
end
% save('BER_MMSE.mat','BER_MMSE');
% 
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR,BER_BDMMSE,'Color','red','LineStyle','-','Marker','^');
xlabel('SNR(in dB)');
ylabel('BER');
legend('MMSEdetect','BDMMSEdetect');