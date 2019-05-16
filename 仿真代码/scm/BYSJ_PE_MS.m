% run('BYSJ_channelmodel.m');
%%
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
for k = 1:7
    SNR(k) = -20+k*2;
    sigma2 = 10^(SNR(k)/10);
    %加噪
    signal_r =awgn(signal_c,SNR(k));
    %%
    % PE检测

    signal_PE = zeros(NoSamples,Nu*Nt);
    H_k =sum(H,3);
    for n_sample=1:NoSamples
            H_A = H_k(:,:,1,n_sample,1);
            for n_link=2:Nu
                H_A =[H_A;H_k(:,:,1,n_sample,n_link)];
            end
                B=H_A*H_A';
                mu =zeros(2*L,1);
                for u=1:2*L
                    mu(u) = trace(B^u)/Nr;
                end
                phi=zeros(L,L);
                for p =1:L
                    for q=1:L
                        phi(p,q) = mu(p+q) + mu(p+q-1)/sigma2;
                    end
                end
                b = (phi)\mu(1:L);
                W=zeros(Nr,Nr);
                for i=1:L
                    W = W + b(i)*(H_A'*H_A)^(i-1);
                end
                W = W *H_A';
                signal_PE(n_sample,:) = W' *signal_r(:,n_sample);
    end
    
    
    %多流判决
     rece_PE = zeros(NoSamples,Nu,Nt);
     for n_link = 1:Nu
         for n_sample = 1:NoSamples
              for n_trans = 1:Nt
                  rece_PE(n_sample,n_link,n_trans) = signal_PE(n_sample,Nt*(n_link-1)+n_trans);
              end
          end
      end
    result_PE = ones(NoSamples,Nu,Nt);
    for i = 1:NoSamples
        for j = 1:Nu
            for p = 1:Nt
                d_1 = abs(rece_PE(i,j,p)-1);
                d_2 = abs(rece_PE(i,j,p)+1);
                if(d_1 > d_2)
                    result_PE(i,j,p) = 0;
                end
            end
        end
    end
    errortimes = sum(abs(result_PE-sample),'all');
    BER_PE(k) = sum(errortimes)/(Nu*Nt*NoSamples);
    
    
end
% save('BER_MMSE.mat','BER_MMSE');
% 
% semilogy(SNR,BER_PE,'Color','blue','LineStyle','-','Marker','o');
% xlabel('SNR');
% ylabel('BER');
% legend('PEdetect');