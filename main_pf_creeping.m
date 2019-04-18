% close all
% clc
% clear all
% April 15 2019
% Yanbing Wang
% Vanderbilt University

%% parameters setup
% model parameters
true = model_param; true.model_name = 'True model';
creep = model_param; creep.model_name = 'Creeping model'; % use model name as key
% lwr = model_param; lwr.model_name = 'LWR';
model_true = true;
model_est = creep;
% vm1 = [1.7*0.9, 1.7*0.95, 1.7, 1.7*1.05, 1.7*1.1];
% rm1 = [1.4*0.9, 1.4*0.95, 1.4, 1.4*1.05, 1.4*1.1];
% rm2 = [1.0*0.9, 1.0*0.95, 1.0, 1.0*1.05, 1.0*1.1];
vm1 = [1.9*0.9, 1.9*0.95, 1.9, 1.9*1.05, 1.9*1.1];
rm1 = [1.6*0.9, 1.6*0.95, 1.6, 1.6*1.05, 1.6*1.1];
rm2 = [1.0*0.9, 1.0*0.95, 1.0, 1.0*1.05, 1.0*1.1];
test = 1; % 1: overtaking, 2: creeping

directory = pwd;
foldername = 'forward_sim';
mkdir(foldername);
directory = fullfile(directory,foldername);

%% initial conditions and boundary conditions
U0_true = initialize(model_true, test);
U0_est = initialize(model_est, test);
% for p1 = 1:5 
% model_est.vm1 = vm1(p1);
% fprintf('vm1 = %.2f',model_est.vm1);
% PF parameters
pf = params_PF;
pf.meas_pt = [3 model_est.N/2 model_est.N-3];
pf.Np = 1500;
pf.init_stdev = 0.06;
pf.model_stdev = 0.10;
pf.meas_stdev = 0.08;
U_true_c1 = zeros([model_true.M model_true.N]); U_true_c2 = U_true_c1;
U_est_c1 = U_true_c1; U_est_c2 = U_true_c1; U_sim_c1 = U_true_c1; U_sim_c2 = U_true_c1; 
U_true_c1(1,:) = U0_true(1,:); U_true_c2(1,:) = U0_true(2,:);
U_est_c1(1,:) = U0_est(1,:); U_est_c2(1,:) = U0_est(2,:);

%% simulation of the true state and the estimator
U_true = cell([1 model_true.M]); U_est = cell([1 model_true.M]); 
U_meas_true = cell([1 model_true.M]); U_res = cell([1 model_true.M]);
U_true{1} = U0_true; U_est{1} = U0_est;
U_meas_true{1} = measure_true(U_true{1},pf);

for n = 1:model_true.M-1
    U_true{n+1} = solver_flow(n,U_true{n},model_true); 
    U_true_c1(n+1,:) = U_true{n+1}(1,:);
    U_true_c2(n+1,:) = U_true{n+1}(2,:);
    
    U_est{n+1} = solver_flow(n,U_est{n},model_est); 
    U_sim_c1(n+1,:) = U_est{n+1}(1,:);
    U_sim_c2(n+1,:) = U_est{n+1}(2,:);
    U_meas_true{n+1} = measure_true(U_true{n+1},pf);
%     if mod(n,10)==0
%         fig = plot_compare(n,U_true,U_est,model_est);
% %         filename = sprintf('simulation_%d_%03d',test,n);
% %         path = fullfile(directory,filename);
% %         saveas(gca,path,'png')
%     end
%     drawnow
end

%% initialize particles and weights
R = covariance(pf,model_est.N);
tau = 0:1:model_est.N-1;
tau_m = triu(toeplitz(tau));
tau_m = tau_m-tau_m';
R_m = autocorrelate(tau_m,model_est.N*50);
[V,D] = eig(R_m);
x_miu = zeros(size(tau));
y_miu = zeros(size(measure(U_est{1},pf)));
init_miu = zeros(size(x_miu));
x = zeros([size(U0_est),pf.Np]);
wt = ones(pf.Np, 1)/pf.Np;
for p = 1:pf.Np
    sum_init = 0;
    for i = 1:size(V,1)
        sum_init = sum_init + sqrt(D(i,i)) * randn * pf.init_stdev * V(:,i)';
    end
    noise_init = init_miu + sum_init; % a random field at time n
    noise_init(noise_init<0) = 0;
    x(:,:,p) = U0_est + noise_init;
    x(x<0)=0;
    %     x(p,1) = rho_est(1,1) + normrnd(0,pf.bound_stdev(1)); % upstream boundary
    %     x(p,end) = rho_est(1,end) + normrnd(0,pf.bound_stdev(2)); % downstream boundary
end
clear sum;

%% Perform time propagation to obtain priori particles
tic
x_next = x; y_next = zeros(size(measure(U_est{1},pf))); N_eff = zeros(1,model_true.M);
U_res{1} = U0_est;
%   ****************************** PARTICLE FILTER STARTS ***********************************
for n = 2:model_true.M
    fprintf('running time = %d\n',n);
    x(x<0)=0;

    d1l_temp = model_est.d1l(n);
    d2l_temp = model_est.d2l(n);
%     d1r_temp = model_est.d1r(n);
%     d2r_temp = model_est.d2r(n);
    for p = 1:pf.Np
        sum_x = 0; sum_y = 0;
        for i = 1:size(V,1)
            sum_x = sum_x + sqrt(D(i,i)) * randn * pf.model_stdev * V(:,i)'; 
        end
        noise_x = x_miu + sum_x;
        noise_y = y_miu + normrnd(0,pf.meas_stdev,size(y_miu));
%------ model predict -----------
        model_est.d1l(n) = d1l_temp + wblrnd(0.06,4);
        model_est.d2l(n) = d2l_temp + wblrnd(0.06,4);
%         model_est.d1r(n) = d1r_temp + wblrnd(0.08,5);
%         model_est.d2r(n) = d2r_temp + wblrnd(0.08,5);
        x_next(:,:,p) = solver_flow(n-1,x(:,:,p),model_est) + noise_x;
        x_next(x_next<0)=0;
%------ measurement update ---------
        y_next(:,:,p) = measure(x_next(:,:,p),pf) + noise_y;
        y_next(y_next<0) = 0;
        y_true = [U_meas_true{n}(1,:)';U_meas_true{n}(2,:)']; % true measurement
        h = [y_next(1,:,p)'; y_next(2,:,p)']; % measurement equation
        m = size(h,1)*size(h,2);
        wt(p) = (2 * pi)^(-m/2) * (sqrt(sum(sum(abs(R).^2))))^(-1/2) * ...
            exp(-1/2 * (y_true - h)'* R^(-1) * (y_true - h));
    end
    
%--- normalize weight ----
    wt = wt./sum(wt); % make sure all the weights of each state i sum up to one
    N_eff(n) = (sum(wt.^2))^(-1); % effective particle size [HOW TO AVOID THE CURSE OF DIMENSIONALITY: SCALABILITY OF PARTICLE FILTERS WITH AND WITHOUT IMPORTANCE WEIGHTS?]
%     disp(N_eff(n))
    
%--- resampling ----
    for p = 1:pf.Np
        x_next(:,:,p) = x_next(:,:,(find(rand <= cumsum(wt),1)));
    end
    %   ******************************* PARTICLE FILTER ENDS **********************************
            
    U_est_temp = mean(x_next,3);
    U_est_temp(U_est_temp<0) = 0;
    U_res{n} = U_est_temp;
    U_est_c1(n,:) = U_est_temp(1,:);
    U_est_c2(n,:) = U_est_temp(2,:);

%     --------------------------- plot and save -------------------
    if mod(n,10)==0
   
        fig = plot_est(n,U_true,U_res,model_est,pf,U_meas_true,x_next);
        drawnow
%     filename = sprintf('pf_test%d_%03d',test,n);
%     path = fullfile(directory,filename);
%     saveas(gca,path,'png')
    end
    x = x_next;
end
toc

%% write video
% write_video(directory,foldername);

%% plot estimated densities
% plot_six(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true);

%% plot effective particle size
% figure; plot(N_eff(2:end),'LineWidth',2,'color',[0.8,0.61,0]); 
% xlabel('Simulation time step','interpreter', 'latex'); ylabel('Effective Particle Size','interpreter', 'latex');
% set(gca,'FontSize',20)
% a = get(gca,'XTick');  
% set(gca,'XTick',a,'FontName','Times');

%% compute MAE
% between simulated and true
MAE_sim_c1 = sum(sum(abs(U_true_c1 - U_sim_c1)))/(model_true.N*model_true.M);
MAE_sim_c2 = sum(sum(abs(U_true_c2 - U_sim_c2)))/(model_true.N*model_true.M);
% between estimated and true
MAE_est_c1 = sum(sum(abs(U_true_c1 - U_est_c1)))/(model_true.N*model_true.M);
MAE_est_c2 = sum(sum(abs(U_true_c2 - U_est_c2)))/(model_true.N*model_true.M);
MAE = [MAE_sim_c1;MAE_sim_c2;MAE_est_c1;MAE_est_c2]
% end