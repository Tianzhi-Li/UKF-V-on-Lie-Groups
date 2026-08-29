% --------------------------------------------------
% Variational Unscented Kalman Filter on SO(3)
% Ref: [Variational Unscented Kalman Filter on Matrix Lie Groups, Tianzhi Li and Jinzhi Wang, Automatica, 172:111995, 2025]
% --------------------------------------------------
function Variational_UKF_SO3()
M1 = 50;     
K1 = M1*(0.1*0.1 + 0.3405*0.3405)/12;   
K2 = M1*(0.1*0.1 + 0.3405*0.3405)/12;    
K3 = M1*(0.1*0.1)/6;   
K = blkdiag(K1,K2,K3);
frequency = 100;        
delta_t = 1/frequency;     
Q = 0.01*eye(3);
N = (pi/180)*(pi/180)*eye(3);
T_long = 50;           
NN = T_long*frequency; 

function output2 = slv_Dynamics(dt,xi,v_xi_post)
    iter = 0;
    maxT = 1e3;         
    tole = 1e-9;   
    x = v_xi_post;       
    err = 1;            
    Jm = zeros(3,3);          
    while err > tole
        Jm(1:3,1:3) = [K1/dt, -0.5*( K2 - K3 )*x(3), -0.5*( K2 - K3 )*x(2);
                       -0.5*( K3 - K1 )*x(3), K2/dt, -0.5*( K3 - K1 )*x(1);
                       -0.5*( K1 - K2 )*x(2), -0.5*( K1 - K2 )*x(1), K3/dt ];       
        if cond(Jm) >= 1e4
            warning('%s','The condition number of the Jacobian matrix is too large')
        end
        xstep = Jm \ Dynamics(dt,xi,x);
        x = x - xstep;
        if iter > maxT
            warning('%s','Exceed maximum iteration steps')
        end
        err = norm( Dynamics(dt,xi,x) );
        iter = iter + 1;
    end
    output2 = x;
end  

function output1 = Dynamics(dt,xi,xi1)                 
    output1 = (1/dt) * K * ( xi1(1:3) - xi(1:3) ) ...
              - 0.5 * Bracket_Star_so3( xi1(1:3) , K*xi1(1:3) ) ...
              - 0.5 * Bracket_Star_so3( xi(1:3) + xi(4:6) , K*(xi(1:3) + xi(4:6) ) );
end

function [v_xi_pri, P_pri] = TimeUpdate(dt,P_post,v_xi_post)
    s1 = 1e-3;
    s2 = 2;
    m = 6;
    lambda_1 = (s1*s1-1)*m;
    w_beta_1 = zeros(13,1);
    w_beta_1(1,1) = lambda_1/( lambda_1+m );
    w_gamma_1 = zeros(13,1);
    w_gamma_1(1,1) = lambda_1/( lambda_1+m ) + 1 - s1*s1 + s2;
    for jj = 2:13                                       
        w_beta_1(jj,1) = 0.5/( lambda_1+m );
        w_gamma_1(jj,1) = 0.5/( lambda_1+m );
    end
    mu_post = zeros(6,13);
    mu_pri = zeros(6,13);
    mu_mean_pri = zeros(3,1);
    P_pri = zeros(3,3);
    P_a_post = blkdiag( P_post, Q );
    Sigma_post = ( chol( (m+lambda_1)*P_a_post ) ).';       
    mu_post(:,1) = zeros(6,1);                        
    for iii = 2:7                                   
        mu_post(:,iii) = Sigma_post(:,iii-1);
        mu_post(:,iii+6) = -Sigma_post(:,iii-1);
    end
    for jjj = 1:13
        mu_pri(1:3,jjj) = slv_Dynamics(dt,mu_post(:,jjj),v_xi_post);
    end
    for jjjj = 1:13  
        mu_mean_pri = mu_mean_pri + w_beta_1(jjjj,1)*mu_pri(1:3,jjjj);   
    end
    for tt = 1:13     
        P_pri = P_pri + w_gamma_1(tt,1)*( mu_pri(1:3,tt) - mu_mean_pri )*(( mu_pri(1:3,tt) - mu_mean_pri ).');   
    end
    v_xi_pri = v_xi_post + mu_mean_pri;
end

function [v_xi_post_post, P_post_post] = MeasurementUpdate(v_xi_pri_2,P_pri_2,z_measure)
    r = 6;
    s1 = 1e-3;
    s2 = 2;
    lambda_2 = (s1*s1-1)*r;
    w_beta_2 = zeros(13,1);
    w_beta_2(1,1) = lambda_2/( lambda_2+r );
    w_gamma_2 = zeros(13,1);
    w_gamma_2(1,1) = lambda_2/( lambda_2+r ) + 1 - s1*s1 + s2;
    for j_j = 2:13                                       
        w_beta_2(j_j,1) = 0.5/( lambda_2+r );
        w_gamma_2(j_j,1) = 0.5/( lambda_2+r );
    end
    eta_pri = zeros(6,13);
    P_xz = zeros(3,3);
    P_zz = zeros(3,3);
    zeta = zeros(3,13);
    zeta_mean = zeros(3,1);
    P_a_pri = blkdiag( P_pri_2, N );
    Sigma_pri = ( chol( (r+lambda_2)*P_a_pri ) ).';        
    eta_pri(:,1) = zeros(6,1);                        
    for iii = 2:7                                      
        eta_pri(:,iii) = Sigma_pri(:,iii-1);
        eta_pri(:,iii+6) = -Sigma_pri(:,iii-1);
    end
    for i_i = 1:13
        zeta(1:3,i_i) = eta_pri(1:3,i_i) + eta_pri(4:6,i_i);     
    end
    for ij = 1:13       
        zeta_mean = zeta_mean + w_beta_2(ij,1)*zeta(1:3,ij);     
    end
    for ii_i = 1:13     
        P_xz = P_xz + w_gamma_2(ii_i,1)*eta_pri(1:3,ii_i)*( ( zeta(:,ii_i) - zeta_mean ).' );    
        P_zz = P_zz + w_gamma_2(ii_i,1)*( zeta(:,ii_i) - zeta_mean )*( ( zeta(:,ii_i) - zeta_mean ).' );  
    end
    K_Gain = P_xz/P_zz;
    P_post_post = P_pri_2 - K_Gain*P_zz*(K_Gain.');
    v_xi_post_post = v_xi_pri_2 + K_Gain*( z_measure - v_xi_pri_2 );  
end

function [outSO3] = so3_Cayley(lambda,dt,xit111)   
    % hat_xit is a 3*3 matrix，xit is a 3D vector
    hat_xit111 = [0,      -xit111(3),  xit111(2);
                  xit111(3),     0,    -xit111(1);
                  - xit111(2),  xit111(1),     0 ];
    cayXi = ( eye(3) - 0.5*dt*hat_xit111 ) \ ( eye(3) + 0.5*dt*hat_xit111 );
    g = lambda;
    output= g * cayXi;        
    outSO3= output(1:3,1:3);        
end



% --------------------------------------------------------------
% Main 
% --------------------------------------------------------------
Reference_omega = zeros(3,NN);
Reference_velocity = zeros(3,NN);
for ii = 1:NN
    Reference_omega(:,ii) = [0.15;0; 0];
    Reference_velocity(:,ii) = Reference_omega(:,ii);
end

Measure_noised = zeros(3,NN);
N_sqrt = chol(N);
Noise = N_sqrt*randn(3,NN);
for iiiii = 1:NN
    Measure_noised(:,iiiii) = Reference_velocity(:,iiiii) + Noise(:,iiiii);
end

Measure_noised_hat = zeros(3,3,NN);
SO3_noised = zeros(3,3,NN);
SO3_noised(:,:,1) = eye(3);
for i1 = 2:NN
    Measure_noised_hat(:,:,i1) = [0, -Measure_noised(3,i1), Measure_noised(2,i1);
                                  Measure_noised(3,i1), 0, -Measure_noised(1,i1);
                                  -Measure_noised(2,i1), Measure_noised(1,i1), 0 ];
    SO3_noised(:,:,i1) = SO3_noised(:,:,i1-1)*expm( delta_t*Measure_noised_hat(:,:,i1-1) );
end

P_priori = zeros(3,3,NN);
v_xi_priori = zeros(3,NN);
P_posterior = zeros(3,3,NN);
P_posterior(:,:,1) = 0.01*eye(3);
v_xi_posterior = zeros(3,NN);  
v_xi_posterior(:,1) = zeros(3,1);    
for kk = 1:(NN-1)    
    [ v_xi_priori(:,kk), P_priori(:,:,kk) ] = TimeUpdate(delta_t,P_posterior(:,:,kk),v_xi_posterior(:,kk) );
    [ v_xi_posterior(:,kk+1), P_posterior(:,:,kk+1) ] = MeasurementUpdate(v_xi_priori(:,kk),P_priori(:,:,kk),Measure_noised(:,kk+1) );
end

SO3 = zeros(3,3,NN);              
SO3(:,:,1) = eye(3);                           
for ii = 1:NN-1        
	SO3(:,:,ii+1) = so3_Cayley(SO3(:,:,ii),delta_t, v_xi_posterior(:,ii+1));
end

err_velocity_omega = zeros(NN,1);

for kkk = 1:NN
    err_velocity_omega(kkk,1) = norm(v_xi_posterior(1:3,kkk) - Reference_velocity(1:3,kkk));      
end

figure
ss = 1:NN;
semilogy( ss*delta_t, err_velocity_omega(ss,1), 'm-', 'LineWidth', 2);
hold on;
ylabel( '$\omega_{k+\frac{1}{2}}^{err}$', 'Interpreter', 'latex' )
title( 'Error in Body Angular Velocity', 'Interpreter', 'latex' )
xlabel( 'Time t (s)', 'Interpreter', 'latex' )

% -------------------------------------
% Ref vel, Measured vel, Est vel
% -------------------------------------
figure
ss = 1:NN;
title( 'Measured, Ground Truth, Estimated', 'Interpreter', 'latex' )
subplot(3,1,1)
plot( ss*delta_t, Measure_noised(1,ss), 'b-' );      % Mea Vel
hold on;
plot( ss*delta_t, Reference_velocity(1,ss), 'c-', 'LineWidth', 3 );  % Ref Vel
hold on;
plot( ss*delta_t, v_xi_posterior(1,ss), 'm-', 'LineWidth', 3 );      % Est Vel
hold on;
subplot(3,1,2)
plot( ss*delta_t, Measure_noised(2,ss), 'b-' );      % Mea Vel
hold on;
plot( ss*delta_t, Reference_velocity(2,ss), 'c-', 'LineWidth', 3 );  % Ref Vel
hold on;
plot( ss*delta_t, v_xi_posterior(2,ss), 'm-', 'LineWidth', 3 );       % Est Vel
hold on;
subplot(3,1,3)
plot( ss*delta_t, Measure_noised(3,ss), 'b-' );       % Mea Vel
hold on;
plot( ss*delta_t, Reference_velocity(3,ss), 'c-', 'LineWidth', 3 );     % Ref Vel
hold on;
plot( ss*delta_t, v_xi_posterior(3,ss), 'm-', 'LineWidth',3 );      % Est Vel
hold on;
xlabel( 'Time t (s)', 'Interpreter', 'latex' )
legend( 'Measured', 'Ground Truth', 'Estimated' )
% -------------------------------------
% Reference traj
% -------------------------------------
time_t = linspace(0,T_long,NN+1);
Reference_R = zeros(3,3,NN);
for ii = 1:NN
    Reference_R(:,:,ii) = [1, 0, 0;
                           0, cos(0.15*time_t(ii)), -sin(0.15*time_t(ii));
                           0, sin(0.15*time_t(ii)), cos(0.15*time_t(ii)) ];
end
% -------------------------------------
% Error in Attitude
% -------------------------------------
err_R = zeros(NN,1);
temp_1 = zeros(3,3,NN);
temp_2 = zeros(3,NN);

for k = 1:NN
    temp_1(:,:,k) = logm( (SO3(:,:,k).')*Reference_R(:,:,k) );
    temp_2(:,k) = [ -temp_1(2,3,k); temp_1(1,3,k); -temp_1(1,2,k) ];    
    err_R(k,1) = norm( temp_2(:,k) );
end
figure
ss = 1:NN;
plot( ss*delta_t, err_R(ss,1), 'm-','LineWidth',2);
ylabel( 'Error in Rotation Matrix', 'Interpreter', 'latex' )
xlabel( 'Time t (s)', 'Interpreter', 'latex' )
xlim([0 50])
ylim([0 1])
grid on;
hold on;

% ---------------------------
% Visualization
% ---------------------------
figure
pic_num1 = 1;
gap = 100;  
x0 = [1;   1;  -1;  -1;   1;  -1;  -1;   1 ];
y0 = [-1;  1;   1;  -1;   1;   1;  -1;  -1 ];
z0 = [-1; -1;  -1;  -1;   1;   1;   1;   1 ];
sec = zeros(8,3);       
sec2 = zeros(8,3); 
for ii = 1:gap:NN
	tmp = SO3(:,:,ii)*([x0 y0 z0]+[0,0,6])';
	sec(:,1) = tmp(1,:)';      
	sec(:,2) = tmp(2,:)';      
	sec(:,3) = tmp(3,:)';      
    tmp2 = SO3_noised(:,:,ii)*([x0 y0 z0]+[0,0,6])';
    sec2(:,1) = tmp2(1,:)';      
	sec2(:,2) = tmp2(2,:)';     
	sec2(:,3) = tmp2(3,:)';      
    RRR = 6;   
    center = [0,0,0];  
    [x, y, z] = ellipsoid(center(1),center(2),center(3),RRR,RRR,RRR,100);  
    surf(x, y, z,'LineStyle','none')
    plot_cubic(sec2)
    colormap winter
    material dull
    axis off    
    M1 = getframe(gcf);
    I1 = frame2im(M1);
    [I1,map] = rgb2ind(I1,256);
    if pic_num1 == 1
        imwrite(I1,map,'UKF-V.gif','gif','Loopcount',inf,'DelayTime',0);
    else
        imwrite(I1,map,'UKF-V.gif','gif','WriteMode','append','DelayTime',0);
    end
    pic_num1 = pic_num1 + 1;
end
hold on
Center = zeros(3,NN);
Center_noised = zeros(3,NN);
for ii = 1:NN
    Center(:,ii) = SO3(:,:,ii)*([0 0 0]+[0,0,6])';
    Center_noised(:,ii) = SO3_noised(:,:,ii)*([0 0 0]+[0,0,6])';
end
sss = 1:NN;
plot3(Center(1,sss),Center(2,sss),Center(3,sss),'LineWidth',3.5,'LineStyle','-.','Color','m')
hold on
hold on



end

