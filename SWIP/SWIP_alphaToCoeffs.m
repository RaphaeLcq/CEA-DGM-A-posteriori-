

function [D] = SWIP_alphaToCoeffs(alpha)

  x_1 = cos(alpha*pi/2);
  y_1 = sin(alpha*pi/2);
  if (abs(y_1) < 1e-16) %% M is always non-invertible, which means any D != 1 will do
    D = 0.5;
    disp('yahou ');
  else
    z_1 = (x_1/y_1);
    z_1_pwr2 = (x_1/y_1)^2;
    D_plus = (2*z_1_pwr2 + 1) + 2*z_1*sqrt(z_1_pwr2 + 1);
    D_min = (2*z_1_pwr2 + 1) - 2*z_1*sqrt(z_1_pwr2 + 1);
    D = min(D_min, D_plus);
  endif
  M = zeros(8,8);
  A = zeros(5,2,2);
  B = zeros(5,2,2);
  for n = 0:4
    theta = n*pi/2;
    x_n = cos(alpha*theta);
    y_n = sin(alpha*theta);
    A(n+1,:,:) = [ x_n , y_n ; -y_n , x_n ];
    B(n+1,:,:) = -[ x_n , y_n ; -D*y_n , D*x_n ];
  endfor
  M(1:2,1:2) = squeeze(A(2,:,:));
  M(1:2,3:4) = squeeze(B(2,:,:));
  M(3:4,3:4) = squeeze(B(3,:,:));
  M(3:4,5:6) = squeeze(A(3,:,:));
  M(5:6,5:6) = squeeze(A(4,:,:));
  M(5:6,7:8) = squeeze(B(4,:,:));
  M(7:8,1:2) = squeeze(A(1,:,:));
  M(7:8,7:8) = squeeze(B(5,:,:));
  global SWIP_coeffs =  robust_null_combined(M)

  SWIP_coeffs = SWIP_coeffs(:,1);
  global DiffusionConstant = D
  fprintf('Diffusion pour alpha = %i doit valoir D = %i\n',alpha,D)
endfunction
