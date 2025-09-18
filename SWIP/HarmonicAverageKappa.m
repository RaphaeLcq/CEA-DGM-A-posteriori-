

function [avg1,avg2] = HarmonicAverageKappa(edge)
  
  
  global EdgTri RefEdg kappa SWIP
  avg1 = 0; avg2 = 0;
  if (SWIP==1)
    if (RefEdg(edge) > 0) %% Dirichlet boundary 
      tri = EdgTri(edge,1);
      avg1 = kappa(tri);
      avg2 = 0;
    else
      tri1 = EdgTri(edge,1); tri2 = EdgTri(edge,2);
      kappa1 = kappa(tri1); kappa2 = kappa(tri2);
      avg1 = (kappa2)/(kappa1 + kappa2);
      avg2 = (kappa1)/(kappa1 + kappa2);
    endif
  else
    avg1 = 0.5;
    avg2 = 0.5;
  endif
  
endfunction
