for i = 1 : size(wave,1) ; imagesc(reshape(wave(i,:),sqrt(size(wave,2)),sqrt(size(wave,2)))) ; drawnow ; colorbar ; end ;


iis = [] ; jjs = [] ; dis = []  ;
MX = 10 ; 
for di = 1 : 100 ; 
   ii = floor((di-1) / MX) + 1 ;
   jj = mod(di, MX) ;
      if(jj==0)
          jj = MX ;
      end
   iis(di) = ii ; jjs(di) = jj ; dis(di) = di ; 
end;