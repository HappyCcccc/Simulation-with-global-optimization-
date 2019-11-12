clear all; close all;

tic
sys_num = 100;
d = 1;
n = 40000;
f = zeros(1,sys_num);
getrho = zeros(1,sys_num);
Z = zeros(1,sys_num);
mu = ones(1,sys_num);
%sm_mu = 1.0;
sigma = ones(1,sys_num);
sigma_squre = ones(1,sys_num);
sum_Z = zeros(1,sys_num);
n_i = zeros(1,sys_num);
sum_sig = zeros(1,sys_num);
sum_sig_squ = zeros(1,sys_num);
n_i_star = zeros(1,sys_num);
Gvn = zeros(1,sys_num);
Mn = 1.0;
current_sys = randi(sys_num,1);
current_fuct = zeros(1,1);

w1 = zeros(1,sys_num);
w2 = zeros(1,sys_num);
disigma =1;
delta_n = zeros(1,n);

a = 241.125;
 for k = 1:sys_num
        f(1,k) = 10*funct(k/sys_num, d);
 end
    
for i = 1:n
    z = disigma * randn;
  
    current_fuct = f(1,current_sys)+ z;
   % current_fuct = current_sys -1 + z;
    %calculate mu: average value of all evalutions in one system
        % sum_Z sum of function value in one system
        
         Z(1,current_sys) = current_fuct;
         sum_Z(1,current_sys) = sum_Z(1,current_sys) + Z(1,current_sys);
         
         if n_i(1,current_sys)>0
             mu(1,current_sys) = sum_Z(1,current_sys)/n_i(1,current_sys);
         end
     
     [sm_mu, n_i_star]= min(mu);
      %calculate sigma 
     w1(1,current_sys) = w1(1,current_sys)+z;
     w2(1,current_sys) = w2(1,current_sys)+z^2;
      
     if n_i(current_sys) >= 2
         sigma(1,current_sys) = sqrt((w2(1,current_sys) - 1/n_i(1,current_sys)*(w1(1,current_sys).^2))/(n_i(1,current_sys)));
         sigma_squre(1,current_sys) = sigma(1,current_sys)^2;
     end    
     
     if n_i(n_i_star) > 0
         Gvn = sqrt(2*log(i*sys_num))*sigma(1,n_i_star)/sqrt(n_i(n_i_star));
     else Gvn =1;
     end
     
    rho = 0.0;
 for j = 1:sys_num
     getrho(j) = sigma_squre(1,j)/(n_i(1,j)*(mu(1,j) - sm_mu + Gvn))^2;     
      if (getrho(j) > rho)
            rho = getrho(j);
            current_sys = j;
      end     
 end  
  n_i(1,current_sys) = n_i(1,current_sys)+1;
  
   p_i  = 1 - n_i(1,n_i_star)/i;
  b_i = pi/4 * sys_num * sqrt(disigma)/sqrt(a)*(2*log(i*sys_num)/i)^0.25;
  rate(1,i) = p_i/b_i;
  
  delta_n(1,i) = (mu(1,n_i_star) - min(f)) * sqrt(i)/(disigma*sqrt(2*log(i*sys_num)));
  
end

plot(delta_n(1,:),'.');
 ylabel({'$\Delta_n\frac{\sqrt{n}}{\sigma_1\sqrt{2\log(n\nu)}}$'},'Interpreter','latex')
 set(gca,'FontSize',20)
  set(gca,'YLim',[-5 2])
 xlabel('Interation number')
 pause(10)

%  plot(rate(1,:),'.');
%  ylabel({'$\Pi_n/B(f,\nu,n)$'},'Interpreter','latex')
%  xlabel('Interation number')
%  pause(30);
x = 1:sys_num;
y = mu;
err = zeros(1,sys_num);
for k=1:sys_num
    err(1,k) = 10*sigma(1,k)/sqrt(n_i(1,k));
end
 
errorbar(x,y,err,'-s','MarkerSize',3,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
axis([0 101 -5 40])

toc
function y = funct(c,d)
%n = size(c);
sum=0;
%sum = (2*c-1.5)^2-cos(6.0*(2*c - 1.5));
sum = 0.5 + 2*c - cos(12*c -9);
y = sum;
end


