for i = 1 : 9 ; 
    clf ; 
    hold on ; 
    bar(2+(i-1)*0.5,1,'r') ;   
    b0 = bar(wave0(i,:),'b','EdgeColor','b') ;  
    b2 = bar(wave2(i,:),'k','EdgeColor','k') ;  
    b0 = bar(wave0(i,:),'w','EdgeColor','b') ;      
    alpha(get(b0,'children'),0);
    set(b0,'LineWidth',2) ;
    set(gca,'XTick',(0.5:0.5:9.5)) ; 
    set(gca,'XGrid','on') ;
    axis([0 9 0 1]); 
    xl = xlabel('i') ; 
    yl = ylabel('\psi_i^N') ; 
    ti = title(['Step N: ' num2str(i-1)]) ; 
    le = legend('Exact solution','Numerical solution, no antidiffusion', 'Numerical solution, 2 antidiffusive iterations') ;  
    set(gca,'XTickLabel',{'','1','','2','','3','','4','','5','','6','','7','','8','','9','','10'}) ;
    set([xl yl ti le gca],'FontSize',18) ; 
    drawnow ;  
    pause ;
end;

clf ; 
subplot(1,2,1) ; 

hold on ; 
p0 = plot(wave0(1,1:100),'r') ; 
set([p0 ],'LineWidth',2) 
xl = xlabel('i') ;
yl = ylabel('\psi_i^{0}') ; 
ti = title('System at start') ; 
set([gca xl yl ti],'FontSize',16) ;

subplot(1,2,2) ;  

hold on ; 
p0 = plot(wave0(140,1:100)) ; 
p1 = plot(wave2(140,1:100),'k') ;
set([p0 p1],'LineWidth',2) 
xl = xlabel('i') ;
yl = ylabel('\psi_i^{140}') ; 
le = legend('No antidiffusion','2 antidiffusive iterations') ;
ti = title('System after 140 iterations') ; 
set([gca xl yl le ti],'FontSize',16) ;


clf ;
hold on ;
p0 = plot(max(wave0'),'k*') ;
p1 = plot(max(wave1'),'r') ;
p2 = plot(max(wave2'),'b') ;
p3 = plot(max(wave3'),'k') ;
set([p0 p1 p2 p3],'LineWidth',2) ;
le = legend('No antidiffusion','1 iteration', '2 iterations', '3 iterations')

ti = title('Maximum value in simulation') ; xl = xlabel('Timestep N') ; yl = ylabel('max(\psi^N)') ;
set([gca le xl yl ti],'FontSize',16) ;



%for i = 1 : 49 ; clf ;subplot(1,2,1)  ; hold on ; bar(10+(i-1)*0.4,1,'r') ; bar(wave2(i,:),'w','EdgeColor','k') ; bar(wave0(i,:),'w','EdgeColor','b') ; subplot(1,2,2) ; bar(antidif(i*2,:),'w','EdgeColor','r') ; drawnow ; pause; end

