sys_num = 30;
 

R = 100;
 

eps = 0.3;
gama = 0.9;
 

original_fun = ones(1, sys_num);
original_mu = zeros(1, sys_num);
original_sigma = zeros(sys_num, sys_num);
 

h = 0;
sum_h = 0;
avg_h;
min_fun;
 

for i = 1 : R
    %initialize
    for j = 1 : sys_num
        for k = 1 : sys_num
            original_sigma(j,k) = exp(-(j-k)^2/sys_num);
        end
    end
    original_fun(1,:)=mvnrnd(original_mu, original_sigma);
    

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
            

end
 

avg_h = norminv(gama^(1/sys_num)) * sum_h/R




