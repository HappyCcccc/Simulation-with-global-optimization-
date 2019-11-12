sys_num = 30;
%samplefun_num = 5;
 

eps = 0.3;
gama = 0.9;
 

R = 1000;
 

original_fun = ones(1, sys_num);
original_mu = zeros(1, sys_num);
original_sigma = zeros(sys_num, sys_num);
 

 

noise = 0;
noise_mu = 0;
noise_sigma = 1;
 

 

current_fun = ones(1, sys_num);
current_mu = ones(1, sys_num);
current_sigma = zeros(sys_num, sys_num);
 

 

total_iteration = ones(1, sys_num);
total_noise = zeros(1, sys_num);
 

interation_count = zeros(1, R);
lte = 0;
gtg = 0;
sum_mu = 0;
avg_mu = 0;
h = 0;
sum_h = 0;
avg_h = 0;
min_fun = 0;
 

E = zeros(sys_num, sys_num);
  

ro = ones(1, sys_num);
 

 

m = zeros(1, sys_num);
 

G = 1;
M = 1;
 
sum_interation = 0;
for i = 1 : R
    %initialize
    for j = 1 : sys_num
        for k = 1 : sys_num
            original_sigma(j,k) = exp(-(j-k)^2/sys_num);
        end
    end
    original_fun(1,:)=mvnrnd(original_mu, original_sigma);
    G = 1;
    ro = ones(1, sys_num);
    total_iteration = ones(1, sys_num);
    total_noise = zeros(1, sys_num);
    

    %caculate h
    min_fun = original_fun(1);
    h = 0;
    for j = 1 : sys_num
        if original_fun(j) < min_fun
            min_fun = original_fun(j);
        end
    end
    for j = 1 : sys_num
        h = h+1/(original_fun(j)-min_fun+eps)^2;
    end
    sum_h = sum_h+h;
            

 

 

    %iteration
    while G >  eps
 

 

    pmu = errorbar(x,y,err,'-s','MarkerSize',3,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
 

    pmu.XDataSource = 'x';
    pmu.YDataSource = 'y';
    pmu.YPositiveDeltaSource = 'err';
    pmu.YNegativeDeltaSource = 'err';
 

    % hold on
    % pmu1 = errorbar(x,z,err,'-s','MarkerSize',3,...
    %     'MarkerEdgeColor','red','MarkerFaceColor','green');
 

    axis([0 31 min(original_fun)-1 max(original_fun)+1])
    set(gca,'FontSize',15)
 

        % update noise
        noise = normrnd(noise_mu, noise_sigma);
 

 

        % interation
        p = 1;
        for j = 1 : sys_num
            if ro(j) > ro(p)
                p = j;
            end
        end
        total_iteration(p) = total_iteration(p)+1;
        total_noise(p) = total_noise(p)+noise; % noise uncertain
 

 

        % update_E
        for j = 1 : sys_num
            E(j, j) = noise_sigma/total_iteration(j);
        end
 

 

        % update_mu
        for j = 1 : sys_num
            m(j) = original_fun(j)+total_noise(j)/total_iteration(j);
        end
        current_mu = original_mu+transpose(original_sigma*inv(original_sigma+E)*transpose(m-original_mu));
 

 

        % update_sigma
        current_sigma = original_sigma-original_sigma*inv(original_sigma+E)*original_sigma;
 

 

        % update_M&G
 

 

        p = 1;
        for j = 1 : sys_num
            if current_mu(j) < current_mu(p)
                p = j;
            end
        end
        M = current_mu(p);
        G = (norminv(gama^(1/sys_num))^2*current_sigma(p, p))^(1/2);
 

 

        % update ro
        for j = 1 : sys_num
            ro(j) = current_sigma(j, j)/(current_mu(j)-M+G)^2;
        end
 

 

    total_iteration;
    current_mu;
 

    x = 1:sys_num;
    y = current_mu;
    err = zeros(1,sys_num);
    z = original_fun;
        for k=1:sys_num
            %err(1,k) = sigma(1,k)/sqrt(n_i(1,k));
            err(1,k) = sqrt(current_sigma(k,k));
        end
 

    % errorbar(x,z,err,'-s','MarkerSize',3,...
    %     'MarkerEdgeColor','green','MarkerFaceColor','green')
    % axis([0 31 -5 5])
    % hold on
 

    refreshdata(pmu);
    drawnow
    

    % count interation times
    interation_count(i) = interation_count(i)+1;
    end
    

    % check if M-min_fun is less than eps
%     sum_mu = 0;
%     for j = 1 : sys_num
%         sum_mu  = sum_mu+current_mu(j);
%     end
%     avg_mu = sum_mu/sys_num;
%     if avg_mu-M < eps
%         lte = lte+1;
%     end


  sum_interation = sum_interation + interation_count(i);



    if M-min_fun < eps
        lte = lte+1;
    end
    

end
 
N=sum_interation/R
lte
H=sum_h/R