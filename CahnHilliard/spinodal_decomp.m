function spinodal_decomp(D,gamma,x,y,tSpan,rc_spline,u0,perimeter_mask,cavity_mask,options)
% Inputs
%    D: double. Diffusion coefficient. Default is 1e-6
%    gamma: double. Square length of transitional regions between domains.
%       Default is 5
% 
% Parameters
%    'GridSize': int. Edge length for the square grid used for the
%       model. Note generation time increases exponentially with grid size.
%       (Default = 200)
%    'FrameSpacing': int. Number of iterations between captured frames, 
%       only relevant if capture mode is standard. (Default = 10)
%    'CaptureMode': char. Method of video cature. Possible inputs below.
%         'standard' - Constant num of iterations between frames. (Default)
%      'incremental' - Num iterations between frames increases over time
%    'ImgStyle': char. Method of frame generation. Possible inputs below.
%         'binary' - Concentrations are binarized to major species. (Default)
%           'true' - Concentrations are mapped to the colormap.
%    'Colormap': char. Colormap used for visualization. Supports all 
%       default MATLAB colormaps (Default = 'jet')
%    'FileName': char. Name of video (Default = 'spinodal_decomposition')
%       
% Examples
%    Example 1: Generate and record default model.
%       spinodal_decomp();
%    Example 2: Generate and record model with custom constants.
%       spinodal_decomp(20,3);
%    Example 3: Generate and record model with custom constants and capture
%    mode.
%       spinodal_decomp(10,10,'CaptureMode','incremental');
%    Example 4: Generate and record model with custom constants and
%    multiple custom parameters.
%       spinodal_decomp(10,10,...
%                      'CaptureMode','incremental',...
%                      'Colormap','jet',...
%                      'ImgStyle','true',...
%                      'NumIterations',25000);   

arguments
   D double
   gamma double
   x (1,:) {mustBeNumeric}
   y (1,:) {mustBeNumeric}
   tSpan (2,1) {mustBeNumeric}
   rc_spline(:,:)
   u0 (:,:) {mustBeNumeric} = 2*randi(2,[100,100])-3;
   perimeter_mask (:,:) = zeros(100,100);
   %oPerimeter (:,:) = zeros(100,100);
   %iPerimeter(:,:) = zeros(100,100);
   cavity_mask(:,:) = ones(100,100);
   options.ImgStyle char = 'true'
   %options.ImgStyle char = 'binary'
   options.Colormap char = 'jet'
   options.GifFileName char = 'EC_phasechange_rc.gif'
   options.ResultsFileName char = 'cahn_hilliard_results.mat'
   options.RunNewSimulation char = 'true'
end

%spinodal_decomp(1e-6,30,'CaptureMode','incremental','dt',60);
% spinodal_decomp(10,10,...
% 'CaptureMode','incremental',...
% 'Colormap','jet',...
% 'ImgStyle','true',...
% 'NumIterations',2500);



% Random Starting Concentrations (-1 and 1 represent different species)
%load circles1.mat s1
%s1(s1==0) = -1;
%u0 = s1;

%u0 = 2*randi(2,options.GridSize)-3;
%u0 = randi(2,options.GridSize)-1;
%x = linspace(0,1,100);
%y = linspace(0,1,100);
%margin = floor(0.2*options.GridSize/2);

%u0 = [-1*zeros(options.GridSize,options.GridSize/2),-1*zeros(options.GridSize,options.GridSize/2)];
%u0 = [-1*ones(options.GridSize,options.GridSize/2),-1*ones(options.GridSize,options.GridSize/2)];
%u0 = -1*ones(options.GridSize,options.GridSize);
%u0(:,options.GridSize/2-20:options.GridSize/2+20) = 2*randi(2,41,options.GridSize)'-3;
%u0(:,floor(options.GridSize/2)-margin:floor(options.GridSize/2)+margin) = -1+2.*rand(options.GridSize,margin*2+1);
%u0(:,options.GridSize/2-20:options.GridSize/2+20) = rand(options.GridSize,41);
%u0(:,options.GridSize/2-20:options.GridSize/2+20) = 0.*ones(options.GridSize,41);
IC = reshape(u0,[],1);

perimeter = mask2ind(perimeter_mask);

[exterior_x,exterior_y] = find(cavity_mask==0);
exterior_inds = sub2ind([length(y),length(x)],exterior_x,exterior_y);

% Getting sparsity matrix for JPattern
S = getSparsityLib.cahn_hilliard_2D(length(x),length(y));
S = sparse(S);

% Setting solver options
opts1 = odeset('JPattern',S,'Vectorized','on','RelTol',1e-3,'AbsTol',1e-4);

% Defining domain (excluding boundaries)
%domain_x = x(2:length(x)-1);
%domain_y = y(2:length(y)-1);
domain_x = 2:length(x)-1;
domain_y = 2:length(y)-1;
domain_x_rep = repmat(domain_x,1,length(domain_y));
domain_y_rep = repelem(domain_y,length(domain_x));
%domain_y_rep = repmat(domain_y,1,length(domain_x));
inner_inds = sub2ind([length(y),length(x)],domain_y_rep,domain_x_rep);
sim_domain_inds = inner_inds;
exterior_common = intersect(sim_domain_inds,exterior_inds);
sim_domain_inds = setxor(sim_domain_inds,exterior_common);
perim_common = intersect(sim_domain_inds,perimeter);
sim_domain_inds = setxor(sim_domain_inds,perim_common);
%[perim_x,perim_y] = ind2sub([length(y),length(x)],perimeter);
[xx,yy] = meshgrid(x,y);
complete_domain = [xx,yy];


% running simulation using ode15s, or just loading the results for plotting
if (options.RunNewSimulation == "true")

    tic
    [t,ut] = ode15s(@(t,u) cahn_hilliard_2D_vec(t,u,x,y,perimeter,sim_domain_inds,D,gamma), tSpan, IC, opts1);
    %[t,ut] = ode15s(@(t,u) cahn_hilliard_2D_vec_rcgrowth(t,u,x,y,complete_domain,inner_inds,rc_spline,D,gamma), tSpan, IC, opts1);
    toc
    
    u = zeros(length(x),length(y),length(t));
    for i = 1:length(t)
        u(:,:,i) = reshape(ut(i,1:length(x)*length(y)),length(x),length(y));
    end

    save(options.ResultsFileName,'t','x','y','u','perimeter','sim_domain_inds')

else
    
    load(options.ResultsFileName)

end

% plotting results and making gif
plotHeatmapGif(t,x,y,u,double(perimeter_mask),'ShowInteriorBounds','true','GifFileName',options.GifFileName)
    
end
