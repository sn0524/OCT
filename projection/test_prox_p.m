clear;
close all

y = rand(200,1024,200);
lam = 0.5;

tic
x1 = prox_psum3(y, y, lam, 3);
disp(['time 1 is ' num2str(toc)]);

tic
x2 = prox_projsum3(y, y, lam, 3);
disp(['time 2 is ' num2str(toc)]);

disp(['diff is ' num2str(norm(x1(:) - x2(:)))]);
