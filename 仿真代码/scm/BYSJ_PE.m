%run('BYSJ_channelmodel.m')
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
signal_c = zeros(Nr,length_sample,Nu);
for n_link=1:Nu
    for n_path=1:NoPath
        for n_sample=1:NoSamples
            H_i = H(:,:,n_path,n_sample,n_link);
            H_i = H_i(:,:);
            signal_c(:,n_sample,n_link) = signal_c(:,n_sample,n_link) + H_i'* signal_t(n_sample,n_link);
        end
    end
end

L = 3;
SNR =zeros(7,1);
BER_PE =zeros(7,1);
for k = 1:7
    SNR(k) = -20+k*5;
    sigma2 = 10^(SNR(k)/10);
    %加噪
    signal_r =awgn(signal_c,SNR(k));
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
                signal_PE(n_sample,n_link) = signal_PE(n_sample,n_link) + W'* signal_r(:,n_sample,n_link);
            end
        end
    end

    %%
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
    BER_PE(k) = sum(errortimes)/(Nu*1000);
end

semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
hold on;
semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
xlabel('SNR');
ylabel('BER');
legend('PEdetect','MMSEdetect');
