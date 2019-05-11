%run('BYSJ_channelmodel.m')
%%
% �ź�����
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
% �����ź�
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


SNR =zeros(7,1);
BER_PE =zeros(7,1);

L = 1;%PE ���ջ�����
for k = 1:7
    SNR(k) = -20+k*2;
    sigma2 = 10^(SNR(k)/10);
    %����
    signal_r =awgn(signal_c,SNR(k));
    %%
    % PE���

    signal_PE = zeros(length_sample,Nu);
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
                W=zeros(Nr,Nt);
                for i=1:L
                    W = W + b(i)*(H_A'*H_A)^(i-1);
                end
                W = W *H_A';
                signal_PE(n_sample,:) = W' *signal_r(:,n_sample);
    end

    %%
    % �о�
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
end
save('BER_PE_3.mat','BER_PE');

% semilogy(SNR,BER_MMSE,'Color','blue','LineStyle','-','Marker','o');
% hold on;
semilogy(SNR,BER_PE,'Color','red','LineStyle','-','Marker','+');
% hold on;
% semilogy(SNR,BER_LPE,'Color','black','LineStyle','-','Marker','*');
% legend('MMSE detect','L=3 PE detect','L=3 LPE detect');
