% Notes_Harker_Les_logo_v2.m
% A small Matlab logo is generated, diff, and then reconstructed with
% Harker 'g2sTikhonov'
%  Uses simplified functions from Harker and O'Leary
%   dopDiffLocal_simple.m
%   dop_simple.m
%   g2sTikhonov_simple.m

%% Generate test surface from Matlab logo
% scaled to 100, gaussian filtered, and with columns = rows - 1
clc; clear;
pathFigures = cd;
% strFileNamePrefix ='logo_micro_';  % micro=5, small=50, large=500 points in membrane
strFileNamePrefix ='logo_small_';  % micro=5, small=50, large=500 points in membrane
% strFileNamePrefix ='logo_large_';  % micro=5, small=50, large=500 points in membrane

membranePoints=50;
if strFileNamePrefix =='logo_micro_', membranePoints = 5; end;
if strFileNamePrefix =='logo_large_', membranePoints = 500; end;
h1 = figure(1);  clf;
Ztrue = 100*membrane(1,membranePoints);  % change the number of points (2nd parameter) over the range 5 to 500
h = fspecial('gaussian',3,1);
Ztrue = imfilter(Ztrue,h);
[m,n]=size(Ztrue)
Ztrue = 100 * Ztrue/max(max(Ztrue));
mean(mean(Ztrue))

s = surface(Ztrue); climZtrue=[0, m, 0, m, -10, 120]; axis(climZtrue);
s.EdgeColor = 'none';
view(3)

Zx = diff(Ztrue,1,2);
Zy = diff(Ztrue,1,1);
m = min(size(Zx))
n=m-1;
Zx = Zx(1:m, 1:n);
Zy = Zy(1:m, 1:n);
Ztrue = Ztrue(1:m, 1:n);

x=linspace(1, n, n)';
y=linspace(1, m, m)';

%%  Add noise to the gradient field 
if strFileNamePrefix =='logo_micro_',
    sigma = 0;
else, 
    sigma = 0.1;
end;
Ax = ( max( Zx(:)) - min( Zx(:) ) )/2 ;
Ay = ( max( Zy(:)) - min( Zy(:) ) )/2 ;
ZxN = Zx + sigma * Ax * randn(m,n) ;
ZyN = Zy + sigma * Ay * randn(m,n) ;

figure(2); clf;
subplot(121); 
s = surface(ZxN);  axis('tight');
    s.EdgeColor = 'none'; 
    h2 =title(['dZ/dx, \sigma = ',num2str(sigma)],'FontSize',16); 
    xlabel('X'); ylabel('Y');
    set(gca,'FontSize',16); 
view(3)

subplot(122);
s = surface(ZyN); axis('tight');
    s.EdgeColor = 'none';
    h2 =title(['dZ/dy, \sigma = ',num2str(sigma)],'FontSize',16); 
    xlabel('X'); ylabel('Y');
    set(gca,'FontSize',16); 
view(3)
print(gcf,[pathFigures,strFileNamePrefix,'gradients.png'],'-dpng')

%%  Tikhonov Regularization
if strFileNamePrefix =='logo_micro_',
    N = 3;  % no effect on Matlab logo
    lambda = 0.025 ; % strong effect, best is  0.025 or smaller
else,
    N = 5;  % no effect on Matlab logo
    lambda = 0.001 ; % strong effect, best is  0.025 or smaller
end;    
deg = 0 ; % strong effect, best is near 0 (non-negative integers only)
Z0 = 0.1*ones(m,n) ;

%% Expanded diff function for Dx
Dx = dopDiffLocal_simple( x, N, N, 'sparse' ) ;
% full(Dx)

Dy = dopDiffLocal_simple( y, N, N, 'sparse' ) ;
% full(Dy)

%%  Simple g2sTikhonov
tic;
[ Ztik, Res ] = g2sTikhonov_simple( ZxN, ZyN, x, y, N, lambda, deg, Z0, Dx, Dy, strFileNamePrefix ) ;
timeRequired=toc

%% Plot results
% Ztik = Z;
h3 = figure(3);  clf;
subplot(121);
s = surface(Ztik + mean(mean(Ztrue - Ztik)));  axis(climZtrue); axis('tight');
    h2 =title(['reconstructed: t=',num2str(timeRequired),' s, \lambda=',num2str(lambda)],'FontSize',16); 
    xlabel('X'); ylabel('Y');
    set(gca,'FontSize',16); 
s.EdgeColor = 'none';
view(3)

subplot(122);
s = surface(Ztrue);  axis(climZtrue);  axis('tight');
    h2 =title('original surface','FontSize',16); 
    xlabel('X'); ylabel('Y');
    set(gca,'FontSize',16); 
s.EdgeColor = 'none';
view(3)
print(gcf,[pathFigures,strFileNamePrefix,'reconstructed_original.png'],'-dpng')

h3 = figure(4);  clf;
residualsMinusMean = Ztrue - Ztik - mean(mean(Ztrue - Ztik));
s = surface(residualsMinusMean);  axis('tight');
    h2 =title('residuals','FontSize',16); 
    xlabel('X'); ylabel('Y');
    set(gca,'FontSize',16); 
s.EdgeColor = 'none';
view(3)
print(gcf,[pathFigures,strFileNamePrefix,'residuals.png'],'-dpng')

%% Display matrices

if strFileNamePrefix =='logo_micro_',
    disp(['original,  ',num2str(size(Ztrue))]);
    disp(Ztrue);
    disp(['Zx ,  ',num2str(size(Zx))]);
    disp(Zx);
    disp(['Zy ,  ',num2str(size(Zy))]);
    disp(Zy);
    disp(['reconstructed ,  ',num2str(size(Ztik))]);
    disp(Ztik + mean(mean(Ztrue - Ztik)));

end;
