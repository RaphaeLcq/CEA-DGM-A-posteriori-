

function mean = meanNeumannEdge(edgGlo)
  
  global NumEdg LgEdg CoorNeu
  
  
  vertices = NumEdg(edgGlo,:);
  XY = CoorNeu(vertices,:);
  [xyp,wp,lambda,np]=IntEdg_Boo5(XY);
  awp = LgEdg(edgGlo)*wp';
  NeuValues = eval_donneesNeumann(xyp(:,1),xyp(:,2),edgGlo);
  mean = sum(awp.*NeuValues);
  mean /= LgEdg(edgGlo);
  
endfunction
