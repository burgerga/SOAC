

for i = 1 : size(wave2,1) ;
    subplot(1,3,1) ;
    imagesc(reshape(wave2(i,:),100,100)) ; caxis([0 4]) ; colorbar ;
    subplot(1,3,2) ;
    imagesc(reshape(mua(i*2,:),101,100)') ; caxis([-2 2]) ;colorbar ;
    subplot(1,3,3) ; 
    imagesc(reshape(mva(i*2,:),100,101)) ; caxis([-2 2]) ;colorbar ;
    drawnow ; pause ;
end
     
