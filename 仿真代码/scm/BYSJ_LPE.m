% run('BYSJ_channelmodel.m')
%
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


SNR =zeros(7,1);
BER_LPE =zeros(7,1);
L = 3;

for k = 1:7
    SNR(k) = -20+k*2;
    sigma2 = 10^(SNR(k)/10);
    %加噪
    signal_r =awgn(signal_c,SNR(k));
    %%
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
                sum_S = sum_S + corrcoef(H_k(:,:,1,1,n_link)'*S(:,:,j,n_link)*H_k(:,:,1,1,n_link));
                if(Nt>1)
                    S(:,:,m,n_link)= S(:,:,m,n_link) + corrcoef(H_k(:,:,1,1,n_link)*E_B(:,:,j)*H_k(:,:,1,1,n_link)')*S(:,:,m-j,n_link);
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
    signal_LPE = zeros(length_sample,Nu);
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
    % 判决
    result_LPE = ones(length_sample,Nu);
    for i = 1:length_sample
        for j = 1:Nu
            d_1 = abs(signal_LPE(i,j)-1);
            d_2 = abs(signal_LPE(i,j)+1);
            if(d_1 > d_2)
                result_LPE(i,j) = 0;
            end
        end
    end
    errortimes = sum(abs(result_LPE-sample));
    BER_LPE(k) = sum(errortimes)/(Nu*NoSamples);
end

save('BER_LPE_3.mat','BER_LPE');
semilogy(SNR,BER_LPE,'Color','red','LineStyle','-','Marker','+');
xlabel('SNR');
ylabel('BER');
legend('LPEdetect');
