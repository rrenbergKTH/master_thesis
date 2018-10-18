function create_bessel_list(min_value_K, max_value_K, N_K, min_value_alpha, max_value_alpha, N_alpha,lambda_vec, Ec_times_beta_vec, factor_vec)
K_values = linspace(min_value_K, max_value_K, N_K);
alpha_values = linspace(min_value_alpha, max_value_alpha, N_alpha);
file_name_meta = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_meta_data.txt';
file_name_0 = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_0.txt';
file_name_1 = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_1.txt';
file_name_2 = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_2.txt';
file_name_3 = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_3.txt';
file_name_4 = '../MCSimulations/MCSimulations/bessel_list_files/bessel_list_4.txt';
fopen(file_name_meta,'w');
dlmwrite(file_name_meta,[Ec_times_beta_vec],'-append','delimiter',' ','roffset',0,'precision',10);
dlmwrite(file_name_meta,[factor_vec],'-append','delimiter',' ','roffset',0,'precision',10);
dlmwrite(file_name_meta,[lambda_vec],'-append','delimiter',' ','roffset',0,'precision',10);
dlmwrite(file_name_meta,[alpha_values],'-append','delimiter',' ','roffset',0,'precision',10);
dlmwrite(file_name_meta,[K_values],'-append','delimiter',' ','roffset',0,'precision',10);

fopen(file_name_0,'w');
Matrix=zeros(length(Ec_times_beta_vec)*length(factor_vec)*N_K,4);
counter = 1;
for s=1:length(Ec_times_beta_vec)
    Ec_times_beta = Ec_times_beta_vec(s);
    for u=1:length(factor_vec)
        factor = factor_vec(u);
        for i=1:N_K
            K = K_values(i);
            Lt = ceil(Ec_times_beta*K*factor);
            Ec_times_delta_tau = Ec_times_beta/Lt;
            bessel_ratio = besseli(1,Ec_times_delta_tau*K^2)/besseli(0,Ec_times_delta_tau*K^2);
            
            if isnan(bessel_ratio) == 1
                bessel_ratio = 1.0;
            end
            value = log(bessel_ratio);
            if isinf(value)
                disp("Warning! Infinite value!");
            end
            Matrix(counter,:) = [Ec_times_beta factor K value];
            counter = counter +1;
            %dlmwrite(file_name_0,[Ec_times_beta factor K bessel_ratio],'-append','delimiter',' ','roffset',0,'precision',10);
        end
    end
end
dlmwrite(file_name_0,Matrix,'-append','delimiter',' ','roffset',0,'precision',10);
disp('25%')
fopen(file_name_1,'w');
Matrix=zeros(length(Ec_times_beta_vec)*length(factor_vec)*length(lambda_vec)*N_K,5);
counter = 1;
for s=1:length(Ec_times_beta_vec)
    Ec_times_beta = Ec_times_beta_vec(s);
    for u=1:length(factor_vec)
        factor = factor_vec(u);
        for j=1:N_K
            K = K_values(j);
            for l=1:length(lambda_vec)
                lambda = lambda_vec(l);
                Lt = ceil(Ec_times_beta*lambda*K*factor);
                Ec_times_delta_tau = Ec_times_beta/Lt;
                bessel_ratio = besseli(1,Ec_times_delta_tau*K^2)/besseli(0,Ec_times_delta_tau*K^2);
                if isnan(bessel_ratio) == 1
                    bessel_ratio = 1.0;
                end
                value = log(bessel_ratio);
                if isinf(value)
                    disp("Warning! Infinite value!");
                end
                %dlmwrite(file_name_1,[Ec_times_beta factor lambda K bessel_ratio],'-append','delimiter',' ','roffset',0,'precision',10);
                Matrix(counter,:) = [Ec_times_beta factor lambda K value];
                counter = counter +1;
            end
        end
    end
end
dlmwrite(file_name_1,Matrix,'-append','delimiter',' ','roffset',0,'precision',10);
disp('50%')
fopen(file_name_2,'w');
Matrix=zeros(length(Ec_times_beta_vec)*length(factor_vec)*N_alpha*N_K,5);
counter = 1;
for s=1:length(Ec_times_beta_vec)
    Ec_times_beta = Ec_times_beta_vec(s);
    for u=1:length(factor_vec)
        factor = factor_vec(u);
        for i=1:N_alpha
            alpha = alpha_values(i);
            for j=1:N_K
                K = K_values(j);
                Lt = ceil(Ec_times_beta/(2*pi*alpha)*factor);
                Ec_times_delta_tau = Ec_times_beta/Lt;
                bessel_ratio = besseli(1,Ec_times_delta_tau*K^2)/besseli(0,Ec_times_delta_tau*K^2);
                if isnan(bessel_ratio) == 1
                    bessel_ratio = 1.0;
                end
                value = log(bessel_ratio);
                if isinf(value)
                    disp("Warning! Infinite value!");
                end
                %dlmwrite(file_name_2,[Ec_times_beta factor alpha K bessel_ratio],'-append','delimiter',' ','roffset',0,'precision',10);
                Matrix(counter,:) = [Ec_times_beta factor alpha K value];
                counter = counter +1;
            end
        end
    end
end
dlmwrite(file_name_2,Matrix,'-append','delimiter',' ','roffset',0,'precision',10);
disp('75%')
fopen(file_name_3,'w');
Matrix=zeros(length(Ec_times_beta_vec)*length(factor_vec)*N_alpha*N_K*length(lambda_vec),6);
counter = 1;
for s=1:length(Ec_times_beta_vec)
    Ec_times_beta = Ec_times_beta_vec(s);
    for u=1:length(factor_vec)
        factor = factor_vec(u);
        for i=1:N_alpha
            alpha = alpha_values(i);
            for j=1:N_K
                K = K_values(j);
                for l=1:length(lambda_vec)
                    lambda = lambda_vec(l);
                    %Lt = ceil(Ec_times_beta/(2*pi*alpha)*lambda^2);
                    Lt = ceil(Ec_times_beta/(2*pi*alpha)*lambda^2*factor);
                    Ec_times_delta_tau = Ec_times_beta/Lt;
                    bessel_ratio = besseli(1,Ec_times_delta_tau*K^2)/besseli(0,Ec_times_delta_tau*K^2);
                    if isnan(bessel_ratio) == 1
                        bessel_ratio = 1.0;
                    end
                    value = log(bessel_ratio);
                    if isinf(value)==1
                        disp("Warning! Infinite value!"); 
                        disp(Lt)
                    end
                    %dlmwrite(file_name_3,[Ec_times_beta factor lambda alpha K bessel_ratio],'-append','delimiter',' ','roffset',0,'precision',10);
                    Matrix(counter,:) = [Ec_times_beta factor lambda alpha K value];
                    counter = counter +1;
                end
            end
        end
    end
end

dlmwrite(file_name_3,Matrix,'-append','delimiter',' ','roffset',0,'precision',10);
disp('100%')

fopen(file_name_4,'w');
Matrix=zeros(length(Ec_times_beta_vec)*length(factor_vec)*N_K,4);
counter = 1;
for s=1:length(Ec_times_beta_vec)
    Ec_times_beta = Ec_times_beta_vec(s);
    for u=1:length(factor_vec)
        factor = factor_vec(u);
        for i=1:N_K
            K = K_values(i);
            Lt = ceil(Ec_times_beta*10*factor);
            Ec_times_delta_tau = Ec_times_beta/Lt;
            bessel_ratio = besseli(1,Ec_times_delta_tau*K^2)/besseli(0,Ec_times_delta_tau*K^2);
            
            if isnan(bessel_ratio) == 1
                bessel_ratio = 1.0;
            end
            value = log(bessel_ratio);
            if isinf(value)
               disp("Warning! Infinite value!"); 
            end
            Matrix(counter,:) = [Ec_times_beta factor K value];
            counter = counter +1;
            %dlmwrite(file_name_0,[Ec_times_beta factor K bessel_ratio],'-append','delimiter',' ','roffset',0,'precision',10);
        end
    end
end
dlmwrite(file_name_4,Matrix,'-append','delimiter',' ','roffset',0,'precision',10);
disp('125%')