function plotHeatmapVideo(t,x,y,u,bounds,options)
    arguments
       t (1,:) {mustBeNumeric}
       x (1,:) {mustBeNumeric}
       y (1,:) {mustBeNumeric}
       u (:,:,:) {mustBeNumeric}
       bounds (:,:) {mustBeNumeric} = zeros(size(u,1),size(u,2));
       options.FrameSpacing double = 2
       options.CaptureMode char = 'standard'
       options.ImgStyle char = 'true'
       %options.ImgStyle char = 'binary'
       options.Colormap char = 'jet'
       options.FileName char = 'EC_phasechange'
       options.CLim (2,1) = [0,1];
    end

    writer = VideoWriter(options.FileName);
    open(writer);

    f = figure();
    if options.ImgStyle == "binary"
        colormap("pink")
    else
        colormap(options.Colormap)
    end
    ax1 = axes("Parent",f);

    imagesc(ax1,x,y,u(:,:,1));
    colorbar(ax1)
    hold on
    ax2 = axes("Parent",f);
    bounds_plot = imagesc(ax2,x,y,bounds);
    linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
    ax2.Visible = 'off';
    ax2.XTick = [];
    xlim(ax1,[-1,1])
    ylim(ax1,[-1,1])
    xlim(ax2,[-1,1])
    ylim(ax2,[-1,1])
    clim(ax1,options.CLim)
    clim(ax2,[0,1])
    alpha(bounds_plot,double(bounds));
    colormap(ax2,"hot");
    frame = getframe(1);
    writeVideo(writer,frame);
   
    for k = 1:length(t)
        if (strcmp(options.CaptureMode,'incremental'))
           if (count == frameStep)
               if (strcmp(options.ImgStyle,'true'))
                   imagesc(ax1,x,y,u(:,:,k));
                   colorbar(ax1)
                   bounds_plot = imagesc(ax2,x,y,bounds);
                   linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
                   ax2.Visible = 'off';
                   ax2.XTick = [];
                   xlim(ax1,[-1,1])
                   ylim(ax1,[-1,1])
                   xlim(ax2,[-1,1])
                   ylim(ax2,[-1,1])
                   clim(ax1,options.CLim)
                   clim(ax2,[0,1])
                   alpha(bounds_plot,double(bounds));
                   colormap(ax2,"hot");
               else
                   uMod = round((u+1)/2); % Binarizes image
                   imagesc(ax1,x,y,uMod(:,:,k));          
                   colorbar(ax1)
                   bounds_plot = imagesc(ax2,x,y,bounds);
                   linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
                   ax2.Visible = 'off';
                   ax2.XTick = [];
                   xlim(ax1,[-1,1])
                   ylim(ax1,[-1,1])
                   xlim(ax2,[-1,1])
                   ylim(ax2,[-1,1])
                   clim(ax1,options.CLim)
                   clim(ax2,[0,1])
                   alpha(bounds_plot,double(bounds));
                   colormap(ax2,"hot");
               end
               frame = getframe(1);
               writeVideo(writer,frame);
               count = 0;
               frameStep = frameStep + 1;
           end
           count = count+1;
           % Standard video mode
           else
               if (mod(k,options.FrameSpacing) == 0)
                   if (strcmp(options.ImgStyle,'true'))
                       imagesc(ax1,x,y,u(:,:,k));
                       colorbar(ax1)
                       bounds_plot = imagesc(ax2,x,y,bounds);
                       linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
                       ax2.Visible = 'off';
                       ax2.XTick = [];
                       xlim(ax1,[-1,1])
                       ylim(ax1,[-1,1])
                       xlim(ax2,[-1,1])
                       ylim(ax2,[-1,1])
                       clim(ax1,options.CLim)
                       clim(ax2,[0,1])
                       alpha(bounds_plot,double(bounds));
                       colormap(ax2,"hot");
                   else
                       uMod = round((u+1)/2);
                       imagesc(ax1,x,y,uMod(:,:,k));
                       colorbar(ax1)
                       bounds_plot = imagesc(ax2,x,y,bounds);
                       linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
                       ax2.Visible = 'off';
                       ax2.XTick = [];
                       xlim(ax1,[-1,1])
                       ylim(ax1,[-1,1])
                       xlim(ax2,[-1,1])
                       ylim(ax2,[-1,1])
                       clim(ax1,options.CLim)
                       clim(ax2,[0,1])
                       alpha(bounds_plot,double(bounds));
                       colormap(ax2,"hot");
                   end
                   frame = getframe(1);
                   writeVideo(writer,frame);
               end
        end
    end
    
    close(writer);

end