function plotHeatmapGif(t,x,y,u,bounds,options)
    arguments
       t (1,:) {mustBeNumeric}
       x (1,:) {mustBeNumeric}
       y (1,:) {mustBeNumeric}
       u (:,:,:) {mustBeNumeric}
       bounds (:,:) {mustBeNumeric} = zeros(size(u,1),size(u,2));
       options.ImgStyle char = 'true'
       %options.ImgStyle char = 'binary'
       options.Colormap char = 'jet'
       options.GifFileName char = 'EC_phasechange.gif'
       options.CLim (2,1) = [0,1];
       options.ShowInteriorBounds = 'true'
    end
    
    f = figure();
    if options.ImgStyle == "binary"
        colormap("pink")
    else
        colormap(options.Colormap)
    end
    ax1 = axes("Parent",f);

    imagesc(ax1,x,y,u(:,:,1));
    colorbar(ax1)
    xlim(ax1,[-1,1])
    ylim(ax1,[-1,1])
    clim(ax1,options.CLim)
    title(ax1,sprintf("t = 0 mins"),'FontSize',24);
    set(gca,'FontSize',20)
    hold on
    if options.ShowInteriorBounds == "true"
        ax2 = axes("Parent",f);
        bounds_plot = imagesc(ax2,x,y,bounds);
        linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
        ax2.Visible = 'off';
        ax2.XTick = [];
        xlim(ax2,[-1,1])
        ylim(ax2,[-1,1])
        clim(ax2,[0,1])
        alpha(bounds_plot,double(bounds));
        colormap(ax2,"hot");
    end

   
    for k = 1:length(t)
       if (strcmp(options.ImgStyle,'true'))
           imagesc(ax1,x,y,u(:,:,k));
           colorbar(ax1)
           xlim(ax1,[-1,1])
           ylim(ax1,[-1,1])
           clim(ax1,options.CLim)
           title(ax1,sprintf("t = %.2f mins",t(k)/60),'FontSize',24);
           set(gca,'FontSize',20)
           hold on
           if options.ShowInteriorBounds == "true"
               ax2 = axes("Parent",f);
               bounds_plot = imagesc(ax2,x,y,bounds);
               linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
               ax2.Visible = 'off';
               ax2.XTick = [];
               xlim(ax2,[-1,1])
               ylim(ax2,[-1,1])
               clim(ax2,[0,1])
               alpha(bounds_plot,double(bounds));
               colormap(ax2,"hot");
           end
       else
           uMod = round((u+1)/2);
           imagesc(ax1,x,y,uMod(:,:,k));
           colorbar(ax1)
           xlim(ax1,[-1,1])
           ylim(ax1,[-1,1])
           clim(ax1,options.CLim)
           title(ax1,sprintf("t = %.2f mins",t(k)/60),'FontSize',24);
           set(gca,'FontSize',20)
           hold on
           if options.ShowInteriorBounds == "true"
               ax2 = axes("Parent",f);
               bounds_plot = imagesc(ax2,x,y,bounds);
               linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
               ax2.Visible = 'off';
               ax2.XTick = [];
               xlim(ax2,[-1,1])
               ylim(ax2,[-1,1])
               clim(ax2,[0,1])
               alpha(bounds_plot,double(bounds));
               colormap(ax2,"hot");
           end
       end
       frame = getframe(f);
       [A,map] = rgb2ind(frame2im(frame),256);
       if k == 1
            imwrite(A,map,options.GifFileName,"gif","LoopCount",Inf,"DelayTime",0.05);
       else
            imwrite(A,map,options.GifFileName,"gif","WriteMode","append","DelayTime",0.05);
       end

    end
end
