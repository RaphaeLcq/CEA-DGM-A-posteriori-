

function RefTriRegion()
  
  global CoorBary Nbtri RefTri 
  global SquareSWIP
  
  
  if (SquareSWIP == 1) 
    for tri = 1:Nbtri
      xy = CoorBary(tri,:);
      if (xy(1) < 1/2 && xy(2) < 1/2)
        RefTri(tri) = 3;
      elseif (xy(1) > 1/2 && xy(2) < 1/2)
        RefTri(tri) = 4;
      elseif (xy(1) < 1/2 && xy(2) > 1/2)
        RefTri(tri) = 2;
      else
        RefTri(tri) = 1;
      endif
    endfor
  endif
  
endfunction
