close all;clearvars
EXP_RANGE = 20;
NOISE_FACTOR = 0.1;

A = ReadCppMatrixFromFile('Q1Inv/Gbd.txt');
[M,N] = size(A);
Ah = A';
AhA = A'*A;
I = eye(size(AhA));

b = ReadCppMatrixFromFile('Q1/Ez_sct_d.txt');
b = b(:,1);

b_noise_phase = rand(size(b))*2*pi - pi;
b_noise_mags  = rand(size(b))*max(abs(b))*NOISE_FACTOR;

b_noise = b_noise_mags.*exp(1.0j*b_noise_phase);

b = b+b_noise;

if M<=N
    disp('Underdetermined');
    x0 = Ah*((A*Ah)\b);
elseif M == N
    x0 = A\b;
else
    disp('Overdetermined');
    x0 = (AhA)\(Ah*b);
end
r0 = A*x0-b;

xx = [];
yy = [];
ll = [];


for lambda_exp = 0:1:EXP_RANGE
    lambda_big = 10^lambda_exp;
    lambda_sml = 10^(-lambda_exp);
    
    x_big = (AhA+lambda_big*I)\(Ah*b);
    x_sml = (AhA+lambda_sml*I)\(Ah*b);
    
    res_big = A*x_big - b;
    res_sml = A*x_sml - b;
    xx = [log(norm(res_sml)),xx,log(norm(res_big))];
    yy = [log(norm(x_sml)),yy,log(norm(x_big))];
    ll = [lambda_sml,ll,lambda_big];
%     scatter(log(norm(res_big)),log(norm(x_big)),'bx');hold on;
%     scatter(log(norm(res_sml)),log(norm(x_sml)),'bx');
%     text(log(norm(res_big)),log(norm(x_big)),num2str(lambda_big));
%     text(log(norm(res_sml)),log(norm(x_sml)),num2str(lambda_sml));
end
figure();
rightmost = log(norm(b));
topmost = log(norm(x0));

right_bar_x = [rightmost,rightmost];
right_bar_y = [min(yy),max(yy)];

top_bar_x = [min(xx),max(xx)];
top_bar_y = [topmost,topmost];

plot(right_bar_x,right_bar_y,'k-.');
hold on;
plot(top_bar_x,top_bar_y,'k-.');

if N<M
    left_bar_x = [log(norm(r0)),log(norm(r0))];
    left_bar_y = [min(yy),max(yy)];
    plot(left_bar_x,left_bar_y,'k-.');
end
plot(xx,yy);
text(xx,yy,num2str(ll.','%3.1e'));
xlabel('log(|Ax-b|)');
ylabel('log(|x|)');




dx = xx(2:end)-xx(1:end-1);
dy = yy(2:end)-yy(1:end-1);
curve = 0*dx;
for ii = 1:(length(dx)-1)
    v1 = [dx(ii);dy(ii)];
    v2 = [dx(ii+1);dy(ii+1)];
    cross_term = (v1(1)*v2(2)-v1(2)*v2(1))/norm(v1)/norm(v2);
    curve(ii) = cross_term;
end



[max_curve,opt_idx] = max(curve);

scatter(xx(opt_idx),yy(opt_idx),'kx');
axis image;grid on;
lambda_opt = ll(opt_idx);

w_opt = (AhA + lambda_opt*I)\(Ah*b);

disp(['Optimum solution achieved at lambda = ' num2str(lambda_opt)]);

