%run('BYSJ_channelmodel.m');
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
            signal_c(:,n_sample,n_link) = signal_c(:,n_sample,n_link) + H_i'* repmat(signal_t(n_sample,n_link),Nt,1);
        end
    end
end

SNR =zeros(7,1);
BER_MMSE =zeros(7,1);
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
                signal_MMSE(n_sample,n_link) = signal_MMSE(n_sample,n_link) + mean(W'* signal_r(:,n_sample,n_link));
            end
        end
    end

    % %%
    % % 判决
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
end
save('BER_MMSE.mat','BER_MMSE');

% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% xlabel('SNR');
% ylabel('BER');
% legend('MMSEdetect');