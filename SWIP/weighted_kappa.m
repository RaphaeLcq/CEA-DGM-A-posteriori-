
function [weighted, shared_weighted] =  weighted_kappa(edge,tri,shared_tri)
  
  global mu kappa TriEdg
  edge_loc = find(TriEdg(tri,:) == edge);
  if  mu(tri,edge_loc) == 1
    [weighted, shared_weighted] = HarmonicAverageKappa(edge);
  else
    [shared_weighted, weighted] = HarmonicAverageKappa(edge);
  endif
  
endfunction   
