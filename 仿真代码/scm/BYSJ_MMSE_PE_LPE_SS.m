run('BYSJ_channelmodel.m');
%%
% 信号生成
sample = rand(NoSamples,Nu);
sample = round(sample);
signal_t = sample;
signal_t(signal_t(:,:)==0) = -1;

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



% L=3; %PE阶数
SNR =zeros(7,1);
BER_MMSE =zeros(7,1);
BER_PE =zeros(7,1);
BER_LPE =zeros(7,1);
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
    
    %%
    for L = 1:3
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
                    b = pinv(phi)*mu(1:L);
                    W=zeros(Nr,Nr);
                    for i=1:L
                        W = W + b(i)*(H_A'*H_A)^(i-1);
                    end
                    W = W *H_A';
                    signal_PE(n_sample,:) = W' *signal_r(:,n_sample);
        end

        %单流判决
         mean_PE = zeros(NoSamples,Nu);
         for n_link = 1:Nu
             for n_sample = 1:NoSamples
                  for n_trans = 1:Nt
                      mean_PE(n_sample,n_link) = mean_PE(n_sample,n_link) + signal_PE(n_sample,Nt*(n_link-1)+n_trans)/Nt;
                  end
              end
          end
        result_PE = ones(NoSamples,Nu);
        for i = 1:NoSamples
            for j = 1:Nu
                d_1 = abs(mean_PE(i,j)-1);
                d_2 = abs(mean_PE(i,j)+1);
                if(d_1 > d_2)
                    result_PE(i,j) = 0;
                end
            end
        end
        errortimes = sum(abs(result_PE-sample));
        BER_PE(k,L) = sum(errortimes)/(Nu*NoSamples);
   
       %%
        % LPE检测
        % LPE接收机参数
    
        E_B =zeros(Nr,Nr,2*L+1);
        S=zeros(Nt,Nt,2*L+1,Nu);
        E_B(:,:,1) = eye(Nr);
        for n_link = 1:Nu
            S(:,:,1,n_link) = eye(Nt);
        end
        H_k =sum(H,3);
        for m = 2:2*L+1
            for j= 1:m-1
                sum_S = zeros(Nr,Nr);
                for n_link = 1:Nu
                    sum_S = sum_S + (H_k(:,:,1,1,n_link)'*S(:,:,j,n_link)*H_k(:,:,1,1,n_link));
                    if(Nt>1)
                        S(:,:,m,n_link)= S(:,:,m,n_link) + (H_k(:,:,1,1,n_link)*E_B(:,:,j)*H_k(:,:,1,1,n_link)')*S(:,:,m-j,n_link);
                    else
                        S(:,:,m,n_link)= S(:,:,m,n_link) + (H_k(:,:,1,1,n_link)*E_B(:,:,j)*H_k(:,:,1,1,n_link)')*S(:,:,m-j,n_link);
                    end
                end
                E_B(:,:,m) = E_B(:,:,m) + sum_S*E_B(:,:,m-j);
            end
        end
        mu_a =zeros(2*L,1);
        for m =1:2*L
            mu_a(m) = trace(E_B(:,:,m+1))/Nr;
        end

        phi_a=zeros(L,L);
        for p =1:L
            for q=1:L
                phi_a(p,q) = mu_a(p+q) + mu_a(p+q-1)/sigma2;
            end
        end
        b_a = pinv(phi_a)*mu_a(1:L);

        % LPE检测
        signal_LPE = zeros(NoSamples,Nu*Nt);
        for n_sample=1:NoSamples
            H_A = H_k(:,:,1,n_sample,1);
            for n_link=2:Nu
                H_A =[H_A;H_k(:,:,1,n_sample,n_link)];
            end
            K=zeros(Nr,Nr);
            for i=1:L
                K = K + b_a(i)*(H_A'*H_A)^(i-1);
            end
                W = K *H_A';
            signal_LPE(n_sample,:) = W'* signal_r(:,n_sample);
        end

        %%
        %单流判决
         mean_LPE = zeros(NoSamples,Nu);
         for n_link = 1:Nu
             for n_sample = 1:NoSamples
                  for n_trans = 1:Nt
                      mean_LPE(n_sample,n_link) = mean_LPE(n_sample,n_link) + signal_LPE(n_sample,Nt*(n_link-1)+n_trans)/Nt;
                  end
              end
          end
        result_LPE = ones(NoSamples,Nu);
        for i = 1:NoSamples
            for j = 1:Nu
                d_1 = abs(mean_LPE(i,j)-1);
                d_2 = abs(mean_LPE(i,j)+1);
                if(d_1 > d_2)
                    result_LPE(i,j) = 0;
                end
            end
        end
        errortimes = sum(abs(result_LPE-sample));
        BER_LPE(k,L) = sum(errortimes)/(Nu*NoSamples);
     end
end
% save('BER_MMSE.mat','BER_MMSE');
% save('BER_PE_4.mat','BER_PE');
% save('BER_LPE_4.mat','BER_LPE');

run('BYSJ_plot');
% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_PE(:,1),'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_PE(:,2),'Color','red','LineStyle','--','Marker','+');
% hold on;
% semilogy(SNR,BER_PE(:,3),'Color','red','LineStyle',':','Marker','+');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSEdetect','L=1 PEdetect','L=2 PEdetect','L=3 PEdetect');

% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
% semilogy(SNR,BER_LPE(:,1),'Color','black','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,2),'Color','black','LineStyle','--','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE(:,3),'Color','black','LineStyle',':','Marker','+');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSEdetect','L=1 LPEdetect','L=2 LPEdetect','L=3 LPEdetect');