function [ Z, Res ] = g2sTikhonov_simple( Zx, Zy, x, y, N, lambda, deg, Z0, Dx, Dy, strFileNamePrefix )
%
% Purpose : Computes the Global Least Squares reconstruction of a surface
%   from its gradient field.
%
% Use (syntax):
%   Z = g2sTikhonov( Zx, Zy, x, y, N, lambda, deg )
%   Z = g2sTikhonov( Zx, Zy, x, y, N, lambda, deg, Z0 )
%   [Z,Res] = g2sTikhonov( Zx, Zy, x, y, N, lambda, deg )
%   [Z,Res] = g2sTikhonov( Zx, Zy, x, y, N, lambda, deg, Z0 )
%
% Input Parameters :
%   Zx, Zy := Components of the discrete gradient field
%   x, y := support vectors of nodes of the domain of the gradient
%   N := number of points for derivative formulas (default=3)
%   lambda := either lambda (1x1) or [ lam ; mu ] (2x1), the regularization parameter(s)
%   deg := The order of the differential regularization terms (typically
%       0,1, or 2).  Standard form is deg=0.
%   Z0 := (optional) a-priori estimate of the solution surface
%
% Return Parameters :
%   Z := The reconstructed surface
%   Res := Residuals of the LS term and the Regularization term, typically
%       needed for calculating the Regularization parameter (e.g., L-Curve)
%
% Description and algorithms:
%   The algorithm solves the normal equations of the Least Squares cost
%   function with Tikhonov Regularization terms, formulated by matrix algebra:
%   e(Z) = || D_y * Z - Zy ||_F^2 + || Z * Dx' - Zx ||_F^2
%         + lambda^2 * || Dy^(deg) * ( Z - Z0 ) ||_F^2 + mu^2 * || ( Z - Z0 ) * Dx'^(deg) ||_F^2
%   The normal equations are a rank deficient Sylvester equation which is
%   solved by means of Householder reflections and the Bartels-Stewart
%   algorithm.
%
% References :
%    @inproceedings{
%    Harker2008c,
%       Author = {Harker, M. and O'Leary, P.},
%       Title = {Least Squares Surface Reconstruction from Measured Gradient Fields},
%       BookTitle = {CVPR 2008},
%       Address= {Anchorage, AK},
%       Publisher = {IEEE},
%       Pages = {1-7},
%          Year = {2008} }
%
%    @inproceedings{
%    harker2011,
%       Author = {Harker, M. and O'Leary, P.},
%       Title = {Least Squares Surface Reconstruction from Gradients:
%           \uppercase{D}irect Algebraic Methods with Spectral, \uppercase{T}ikhonov, and Constrained Regularization},
%       BookTitle = {IEEE CVPR},
%       Address= {Colorado Springs, CO},
%       Publisher = {IEEE},
%       Pages = {2529--2536},
%          Year = {2011} }
%
% Author :  Matthew Harker and Paul O'Leary
% Date :    17. January 2013
% Version : 1.0
%
% (c) 2013 Matthew Harker and Paul O'Leary, 
% Chair of Automation, University of Leoben, Leoben, Austria
% email: office@harkeroleary.org, 
% url: www.harkeroleary.org
%
% History:
%   Date:           Comment:
%   Feb. 14, 2011   Original Version
%
    lam = lambda ;
    mu = lambda ;
%====================================================================
% Generate the Differentiation Matrices and Solve the Sylvester Equation
%========================================================================
%
% Dx = dopDiffLocal( x, N, N, 'sparse' ) ;
% Dy = dopDiffLocal( y, N, N, 'sparse' ) ;
% disp('Dx'); full(Dx)
%
tol = sqrt( eps(1) ) ;
%
  %
    A = [ Dy ; mu * speye(length(y)) ] ;
    B = [ Dx ; lam * speye(length(x)) ] ;
    F = [ Zy ; mu * Z0 ] ; % Degree-0 means Dx^0 = I
    G = [ Zx , lam * Z0 ] ;
%=========================
% lyap is in the Matlab Control Systems Toolbox
%=========================
    Z = lyap( full(A'*A), full(B'*B), -A'*F - G*B ) ;

%=========================
% Compute the Residuals:
%=========================
%
    Res = zeros(2,2) ;
    
    Res(1,1) = norm( Z * Dx' - Zx, 'fro' ) ;
    Res(1,2) = norm( Dy * Z  - Zy, 'fro' ) ;
    
if strFileNamePrefix =='logo_micro_',
    disp(['A (converted from sparse to full),  ',num2str(size(A))]);
    disp(full(A));
    disp(['B (converted from sparse to full),  ',num2str(size(B))]);
    disp(full(B));
    disp(['F ,  ',num2str(size(F))]);
    disp(F);    
    disp(['G ,  ',num2str(size(G))]);
    disp(G);    
    disp(['Z ,  ',num2str(size(Z))]);
    disp(Z);    
end;
   
%========
% END
%======== 

end