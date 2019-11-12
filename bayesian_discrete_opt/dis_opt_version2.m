sys_num = 30;
samplefun_num = 5;
 

eps = 0.3;
gama = 0.9;
 

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
 

E = zeros(sys_num, sys_num);
 

ro = ones(1, sys_num);
 

m = zeros(1, sys_num);
 

G = 1;
M = 1;
 

%initialize
for i = 1 : sys_num
    for j = 1 : sys_num
        original_sigma(i,j) = exp(-(i-j)^2/sys_num);
    end
end
original_fun(1,:)=mvnrnd(original_mu, original_sigma);
 

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
    for i = 1 : sys_num
        if ro(i) > ro(p)
            p = i;
        end
    end
    total_iteration(p) = total_iteration(p)+1;
    total_noise(p) = total_noise(p)+noise; % noise uncertain
    

    % update_E
    for i = 1 : sys_num
        E(i, i) = noise_sigma/total_iteration(i);
    end
    

    % update_mu
    for i = 1 : sys_num
        m(i) = original_fun(i)+total_noise(i)/total_iteration(i);
    end
    current_mu = original_mu+transpose(original_sigma*inv(original_sigma+E)*transpose(m-original_mu));
    

    % update_sigma
    current_sigma = original_sigma-original_sigma*inv(original_sigma+E)*original_sigma;
    

    % update_M&G
   
    
    p = 1;
    for i = 1 : sys_num
        if current_mu(i) < current_mu(p)
            p = i;
        end
    end
    M = current_mu(p);
    G = (norminv(gama^(1/sys_num))^2*current_sigma(p, p))^(1/2);
    

    % update ro
    for i = 1 : sys_num
        ro(i) = current_sigma(i, i)/(current_mu(i)-M+G)^2;
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

end



