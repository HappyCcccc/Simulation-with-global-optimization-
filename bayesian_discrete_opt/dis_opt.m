clear all; close all;

tic
sys_num = 30;
samplefun_num = 5;
d = 1;
eps = 0.0000001;
gma =0.9;
f = zeros(1,sys_num);
getrho = zeros(1,sys_num);
Z = zeros(1,sys_num);
mu = ones(1,sys_num);
%sm_mu = 1.0;
sigma = ones(1,sys_num);
sigma_squ = ones(sys_num,sys_num);  % sigma_squ(i,i) = sigma^2/n_i
sum_Z = zeros(1,sys_num);
n = ones(1,sys_num);
sum_sig = zeros(1,sys_num);
sum_sig_squ = zeros(1,sys_num);
n_i_star = zeros(1,sys_num);
Gn = 1;
Mn = 1.0;
current_sys = randi(sys_num,1);
current_fuct = zeros(1,1);
gau_mu = ones(1,sys_num);
gau_sig = zeros(sys_num, sys_num);
update_mu = ones(1,sys_num);
update_sig = zeros(sys_num, sys_num);


for i = 1:sys_num
    for j = 1:sys_num
        gau_sig(i,j) = exp(-(i-j)^2/sys_num);
    end
end

%for j = 1:samplefun_num

%  for k = 1:sys_num
%       %  f(1,k) = 10*funct(k/sys_num, d);
%       f(1,k) = normrnd(0, 1);
%  end

f(1,:)=mvnrnd(gau_mu, gau_sig);

count = 1;     
current_sys = 1;
%while Gn > eps
for i = 1:2000
    %z = randn;
    
    z = normrnd(0,sigma);
  
   % current_fuct = current_sys -1 + z;
    %calculate mu: average value of all evalutions in one system
        % sum_Z sum of function value in one system
        
        sum_Z(1,current_sys) = sum_Z(1,current_sys) + z(1,current_sys)
        mu(1,current_sys) = f(1,current_sys) + sum_Z/n(1,current_sys);
        
      %calculate sigma   E?i?i) = sigma/n_i;
      for k = 1:sys_num
          l = 1;
          while l == k
              sigma_squ(k,l) = sigma^2/n(1,current_sys);
          end
          l = l + 1;
      end
       
 % calculate mu: u_n = gau_mu - gau_sig(gau_sig + f + noise_mean - gau_mu)
  
   update_mu = gau_mu + transpose(gau_sig*inv(gau_sig + sigma_squ)*transpose((mu(1,current_sys)-gau_mu))) ; 
      
   %calculate covarian matrix R_n = R_0 - R_0(R_0+E)^(-1)R_0 
   
   update_sig = gau_sig - gau_sig*inv(gau_sig + sigma_squ)* gau_sig;
      
     [sm_mu, n_i_star]= min(update_mu);
     
     if n(n_i_star) > 0
        % Gn = sqrt(norminv(gma^(1/sys_num))^2* update_sig(n_i_star));
        Gn =1;
     else Gn =1;
     end
     
      if Gn <= eps
          n(1,current_sys) = n(1,current_sys)+1;
         break
       end
     
    rho = 0.0;
 for j = 1:sys_num
     getrho(j) = update_sig(1,j)/(update_mu(1,j) - sm_mu + Gn)^2;     
      if (getrho(j) > rho)
            rho = getrho(j);
            current_sys = j;
      end     
 end  
  n(1,current_sys) = n(1,current_sys)+1;
  count = count+1;
end
% figure;
%   ylabel({'$\Pi_n/B(f,\nu,n)$'},'Interpreter','latex')
%  
x = 1:sys_num;
y = update_mu;
err = zeros(1,sys_num);
z = f;
for k=1:sys_num
    %err(1,k) = sigma(1,k)/sqrt(n_i(1,k));
    err(1,k) =sigma(1,k);
end
errorbar(x,z,err,'-s','MarkerSize',3,...
    'MarkerEdgeColor','green','MarkerFaceColor','green')
axis([0 31 -10 40])
hold on
errorbar(x,y,err,'-s','MarkerSize',3,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
axis([0 31 -10 40])
set(gca,'FontSize',15)
toc

%end



