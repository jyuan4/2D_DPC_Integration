pathFigures = cd;
%% read fits files (H/V bunny)
HoriFilename = [cd, '/Expt9_downsize_withFilter.fits'];
VertFilename = [cd, '/Expt10_downsize_withFilter.fits'];
%HoriFilename = [cd, '/afterCorr_hori_144.fits'];
%VertFilename = [cd, '/afterCorr_vert_144.fits'];
VertData = fitsread(VertFilename);
HoriData = fitsread(HoriFilename);
%%
strFileNamePrefix ='logo_small_';  % micro=5, small=50, large=500 points in membrane
% strFileNamePrefix ='logo_large_'
;  % micro=5, small=50, large=500 points in membrane
membranePoints=50;
if strFileNamePrefix =='logo_micro_', membranePoints = 5; end;
if strFileNamePrefix =='logo_large_', membranePoints = 500; end;
ZxN = 1.5*HoriData(20:230,20:220);
ZyN = VertData(20:230,20:220);

%% cropping
[m, n] = size(ZxN);
x=linspace(1, n, n)';
y=linspace(1, m, m)';
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
[ Ztik, Res ] = g2sTikhonov_simple( ZxN, ZyN, x, y, N, lambda, deg, Z0, Dx, Dy, strFileNamePrefix);
imshow(Ztik,[-4,4])
fitswrite(Ztik, 'phase153_1.5X.fits')

%%
%%
%%
%% Loop
pathForExpt9 = '/Volumes/data2/APS/APS_July2017/Expt9_foramA_GH_6mm/Integration/';
pathForExpt10 = '/Volumes/data2/APS/APS_July2017/Expt10_foramA_GV_6mm/Integration/';
Expt9Files = dir('/Volumes/data2/APS/APS_July2017/Expt9_foramA_GH_6mm/Integration/Expt9*.*');
Expt9Files(1).name
Expt10Files = dir('/Volumes/data2/APS/APS_July2017/Expt10_foramA_GV_6mm/Integration/Expt10*.*');
Expt10Files(1).name

for k = 1: length(Expt9Files)
    VertData = fitsread([pathForExpt10, Expt10Files(k).name]);
    HoriData = fitsread([pathForExpt9,Expt9Files(k).name]);
    
    strFileNamePrefix ='logo_small_'; 
    membranePoints=50;
    if strFileNamePrefix =='logo_micro_', membranePoints = 5; end;
    if strFileNamePrefix =='logo_large_', membranePoints = 500; end;
    ZxN = 2.0*HoriData(20:230,20:220);
    ZyN = VertData(20:230,20:220);

    [m, n] = size(ZxN);
    x=linspace(1, n, n)';
    y=linspace(1, m, m)';
    
    % Tikhonov Regularization
    if strFileNamePrefix =='logo_micro_',
        N = 3;  % no effect on Matlab logo
        lambda = 0.025 ; % strong effect, best is  0.025 or smaller
    else,
        N = 5;  % no effect on Matlab logo
        lambda = 0.001 ; % strong effect, best is  0.025 or smaller
    end;    
    deg = 0 ; % strong effect, best is near 0 (non-negative integers only)
    Z0 = 0.1*ones(m,n) ;
    
    % Expanded diff function for Dx
    Dx = dopDiffLocal_simple( x, N, N, 'sparse' ) ;
    % full(Dx)
    Dy = dopDiffLocal_simple( y, N, N, 'sparse' ) ;
    % full(Dy)'
    
    %  Simple g2sTikhonov
    tic;
    [ Ztik, Res ] = g2sTikhonov_simple( ZxN, ZyN, x, y, N, lambda, deg, Z0, Dx, Dy, strFileNamePrefix);
    %imshow(Ztik,[-4,4])
    FileName = ['Phase_2.0X_Hori9_', num2str(k), '.fits'];
    fitswrite(Ztik, ['/Volumes/data2/APS/APS_July2017/Phase_2.0X_Hori_Exp9/', FileName])
end



