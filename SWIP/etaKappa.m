

function etaK = etaKappa(edge)
  
  global EdgTri RefEdg kappa
  if (RefEdg(edge) > 0) %%  boundary 
    tri = EdgTri(edge,1);
    etaK = kappa(tri);
  else
    tri1 = EdgTri(edge,1); tri2 = EdgTri(edge,2);
    kappa1 = kappa(tri1); kappa2 = kappa(tri2);
    etaK = (2*kappa1*kappa2)/(kappa1 + kappa2);
    
  endif
endfunction
