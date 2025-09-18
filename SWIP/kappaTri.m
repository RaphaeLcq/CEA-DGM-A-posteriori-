

function kappa = kappaTri(D)
  
  global RefTri Nbtri
  kappa = ones(Nbtri,1);
  Ind = find(RefTri == 2 | RefTri == 4)';
  kappa(Ind) = D;
  
  
endfunction
