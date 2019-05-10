run('BYSJ_channelmodel.m');
%%
% 信号生成
length_sample = NoSamples;
sample = rand(length_sample,Nu);
sample = round(sample);
signal_t = ones(length_sample,Nu);
for i = 1:length_sample
    for j = 1:Nu
        if(sample(i,j) == 0)
            signal_t(i,j) = -1;
        end
    end
end

%%
% 接收信号
signal_c = zeros(Nr,length_sample);
for n_link=1:Nu
    for n_path=1:NoPath
        for n_sample=1:NoSamples
            H_i = H(:,:,n_path,n_sample,n_link);
            H_i = H_i(:,:);
            signal_c(:,n_sample) = signal_c(:,n_sample) + H_i'* repmat(signal_t(n_sample,n_link),Nt,1);
        end
    end
end


L=3; %PE阶数

SNR =zeros(7,1);
BER_MMSE =zeros(7,1);
BER_PE =zeros(7,1);
BER_LPE =zeros(7,1);
for k = 1:7
    SNR(k) = -20+k*5;
    sigma2 = 10^(SNR(k)/10);
    %加噪
    signal_r =awgn(signal_c,SNR(k));
    %%
    % MMSE检测

    signal_MMSE = zeros(length_sample,Nu);
    for n_link=1:Nu
        for n_path=1:NoPath
            for n_sample=1:NoSamples
                H_i = H(:,:,n_path,n_sample,n_link);
                H_i = H_i(:,:);
                W = inv(H_i'*H_i+eye(Nr)/sigma2)*H_i';
                signal_MMSE(n_sample,n_link) = signal_MMSE(n_sample,n_link) + mean(W'* signal_r(:,n_sample));
            end
        end
    end


    % 判决
    result_MMSE = ones(length_sample,Nu);
    for i = 1:length_sample
        for j = 1:Nu
            d_1 = abs(signal_MMSE(i,j)-1);
            d_2 = abs(signal_MMSE(i,j)+1);
            if(d_1 > d_2)
                result_MMSE(i,j) = 0;
            end
        end
    end
    errortimes = sum(abs(result_MMSE-sample));
    BER_MMSE(k) = sum(errortimes)/(Nu*NoSamples);
    
    %%
    % PE检测
    signal_PE = zeros(length_sample,Nu);
    for n_link=1:Nu
        for n_path=1:NoPath
            for n_sample=1:NoSamples
                H_i = H(:,:,n_path,n_sample,n_link);
                H_i = H_i(:,:);
                B=H_i*H_i';
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
                W=zeros(Nr,Nt);
                for i=1:L
                    W = W + b(i)*(H_i'*H_i)^(i-1);
                end
                W = W *H_i';
                signal_PE(n_sample,n_link) = signal_PE(n_sample,n_link) + mean(W'* signal_r(:,n_sample));
            end
        end
    end

    % 判决
    result_PE = ones(length_sample,Nu);
    for i = 1:length_sample
        for j = 1:Nu
            d_1 = abs(signal_PE(i,j)-1);
            d_2 = abs(signal_PE(i,j)+1);
            if(d_1 > d_2)
                result_PE(i,j) = 0;
            end
        end
    end
    errortimes = sum(abs(result_PE-sample));
    BER_PE(k) = sum(errortimes)/(Nu*NoSamples);
    
   %%
    % LPE检测
    %%LPE 确定性等同计算
    E_B =zeros(Nr,Nr,2*L+1);
    S=zeros(Nt,Nt,2*L+1,Nu);
    E_B(:,:,1) = eye(Nr);
    for n_link = 1:Nu
        S(:,:,1,n_link) = eye(Nt);
    end

    for m = 2:2*L+1
        for j= 1:m-1
            sum_S = zeros(Nr,Nr);
            for n_link = 1:Nu
                H_k = H(:,:,1,1,n_link);
                H_k = H_k(:,:);
                sum_S = sum_S + corrcoef(H_k'*S(:,:,j,n_link)*H_k);
                if(Nt>1)
                    S(:,:,m,n_link)= S(:,:,m,n_link) + corrcoef(H_k*E_B(:,:,j)*H_k')*S(:,:,m-j,n_link);
                else
                    S(:,:,m,n_link)= S(:,:,m,n_link) + (H_k*E_B(:,:,j)*H_k')*S(:,:,m-j,n_link);
                end
            end
            E_B(:,:,m) = E_B(:,:,m) + sum_S*E_B(:,:,m-j);
        end
    end
    mu_a =zeros(2*L,1);
    for m =1:2*L
        mu_a(m) = real(trace(E_B(:,:,m+1))/Nr);
    end
    phi_a=zeros(L,L);
    for p =1:L
        for q=1:L
            phi_a(p,q) = mu_a(p+q) + mu_a(p+q-1)/sigma2;
        end
    end
    b_a = pinv(phi_a)*mu_a(1:L);
    
    
    signal_LPE = zeros(length_sample,Nu);
    for n_link=1:Nu
        for n_path=1:NoPath
            for n_sample=1:NoSamples
                H_i = H(:,:,n_path,n_sample,n_link);
                H_i = H_i(:,:);
                for i=1:L
                    W = W + b_a(i)*(H_i'*H_i)^(i-1);
                end
                W = W *H_i';
                signal_LPE(n_sample,n_link) = signal_LPE(n_sample,n_link) + mean(W'* signal_r(:,n_sample));
            end
        end
    end

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
save('BER_MMSE.mat','BER_MMSE');

semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
hold on;
semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
legend('MMSE detect','L=3 PE detect','L=3 LPE detect');