load wave.dat
dim = 50
for i = 1 : size(wave,1) ; 
    ma = reshape(wave(i,:),dim, dim) ; 
    clf ; 
    surf(ma) ; 
    %axis([0 dim 0 dim 0 1]) ;
    drawnow ;
end ;