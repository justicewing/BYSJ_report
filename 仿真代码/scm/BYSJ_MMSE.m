run('BYSJ_channelmodel.m');
%%
% 信号生成
sample = rand(NoSamples,Nu);
sample = round(sample);
signal_t = sample;
signal_t(signal_t(:,:)==0) = -1;
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
            signal_c(:,n_sample) = signal_c(:,n_sample) + H_i'* repmat(signal_t(n_sample,n_link),Nt,1);
        end
    end
end
    

SNR =zeros(7,1);
BER_MMSE =zeros(7,1);
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
    
    
    %单流判决
     mean_MMSE = zeros(NoSamples,Nu);
     for n_link = 1:Nu
         for n_sample = 1:NoSamples
              for n_trans = 1:Nt
                  mean_MMSE(n_sample,n_link) = mean_MMSE(n_sample,n_link) + signal_MMSE(n_sample,Nt*(n_link-1)+n_trans)/Nt;
              end
          end
      end
    result_MMSE = ones(NoSamples,Nu);
    for i = 1:NoSamples
        for j = 1:Nu
            d_1 = abs(mean_MMSE(i,j)-1);
            d_2 = abs(mean_MMSE(i,j)+1);
            if(d_1 > d_2)
                result_MMSE(i,j) = 0;
            end
        end
    end
    errortimes = sum(abs(result_MMSE-sample));
    BER_MMSE(k) = sum(errortimes)/(Nu*NoSamples);
    
    
end
save('BER_MMSE2.mat','BER_MMSE');
% 
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
xlabel('SNR');
ylabel('BER');
legend('MMSEdetect');