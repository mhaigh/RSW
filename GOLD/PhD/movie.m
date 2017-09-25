% MOVIE
% File for plotting and saving a movie

% Load the file to be plotted. This file should be 
load /home/mike/Documents/GulfStream/Code/DATA/2L/BALANCED/UNIFORM/u_nd_CENTER128.mat;
load /home/mike/Documents/GulfStream/Code/DATA/2L/BALANCED/UNIFORM/v_nd_CENTER128.mat;
load /home/mike/Documents/GulfStream/Code/DATA/2L/BALANCED/UNIFORM/eta_nd_CENTER128.mat;
load /home/mike/Documents/GulfStream/Code/DATA/1L/BALANCED/UNIFORM/P_NORTH256.mat;

var = 3;

if var == 1  
    
    xdim=size(u_nd,1);
    ydim=size(u_nd,2);
    Nt=size(u_nd,3);
    
    ulim=max(max(max(abs(u_nd))));
    
    for ii=1:4*Nt
        ii
        jj=mod(ii,Nt)+1;
        figure(1)
        surf(u_nd(:,:,jj),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet);...
            caxis([-ulim ulim]);
            
     
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/u.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/u.gif','WriteMode','append');

        end
    end
end

if var == 2 
    
    xdim=size(v_nd,1);
    ydim=size(v_nd,2);
    Nt=size(v_nd,3);
    
    vlim=max(max(max(abs(v_nd))));
    
    for ii=1:4*Nt
        ii
        jj=mod(ii,Nt);        
        figure(1)
        surf(v_nd(:,:,jj),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet);...
            caxis([-vlim vlim]);
     
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/v.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/v.gif','WriteMode','append');

        end
    end
end

if var == 3  
    
    xdim=size(eta_nd,1);
    ydim=size(eta_nd,2);
    Nt=size(eta_nd,3);
    
    etalim=max(max(max(abs(eta_nd))));
    
    for ii=1:4*Nt
        ii
        jj=mod(ii,Nt)+1; 
        figure(1)
        surf(eta_nd(:,:,jj),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet);...
            caxis([-0.8*etalim 0.8*etalim]),axis([1 ydim 1 xdim]);
     
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/eta.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/eta.gif','DelayTime',0,'WriteMode','append');

        end
    end
end

if var == 4  
    
    xdim=size(P,1);
    ydim=size(P,2);
    Nt=size(P,3);
    
    Plim=max(max(max(abs(P))));
    
    for ii=1:4*Nt
        ii
        jj=mod(ii,Nt)+1; 
        figure(1)
        surf(P(:,:,jj),'edgecolor','none'); shading interp; view(0,90); colorbar; colormap(jet);...
            caxis([-Plim Plim]),axis([1 ydim 1 xdim]);
     
        mov(ii)=getframe(1);
        im=frame2im(mov(ii));
        [imind,cm] = rgb2ind(im,256);

        if ii == 1

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/eta.gif', 'Loopcount',inf);

        else

            imwrite(imind,cm,'/home/mike/Documents/GulfStream/Code/IMAGES/eta.gif','DelayTime',0,'WriteMode','append');

        end
    end
end
