classdef getPlotsLib
%     show_strain = 1;
%     show_cavity_radius = 0;
%     show_pressure_field = 1;
%     show_pc_pcrit_p0 = 0;
%     show_p_inf = 0;
%     show_heatmaps = 0;
%     show_concentration_lineplots = 0;
%     show_mass_conservation = 0;

    methods(Static)
        function plot_cavity_radius_cell(t, cavity_radius, param_vec, varargin)
            
            % setting optional parameter values, including parsing for
            % validation. parameters are: 
            % x_lim_end - setting the x_lim end value
            % x_units - setting the units for the x-axis
            default_x_lim_end = t(end);
            default_x_units = 'mins';
            expected_x_units = {'secs','mins','hrs','days'};

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            addParameter(p, 'x_units', default_x_units, ...
                @(x) any(validatestring(x,expected_x_units)));
            
            parse(p,varargin{:});
            
            switch p.Results.x_units
                case 'mins'
                    divisor = 60;
                case 'hrs'
                    divisor = 60*60;
                case 'days'
                    divisor = 24*60*60;
                otherwise
                    divisor = 1;
            end
            
            %plotting
            arrayLegend=cell(length(param_vec),1);
            for i = 1:length(param_vec)
                arrayLegend{i} = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
            end
            
            f = figure();
            ax = axes('Parent', f);

            cm = colormap('copper');
            
            [numColors,~] = size(cm);
            cm_divisor = floor(numColors/length(param_vec));
            
            %[~, col] = size(cavity_radius);
            
            for i = 1:length(cavity_radius)
                timeVec = t{i};
                cavity_radiusVec = cavity_radius{i};
                plot(ax,timeVec/divisor, cavity_radiusVec, 'LineWidth', 4, 'color', cm(cm_divisor*i,:));
                hold on
            end
            x_axis_label = sprintf('Time (%s)', p.Results.x_units);
            xlabel(ax,x_axis_label, 'FontSize', 32)
            ylabel(ax,'Cavity Radius (cm)','FontSize', 32)
            set(ax,'FontSize',28)
            legend(ax,arrayLegend)
            xlim(ax,[0, p.Results.x_lim_end/divisor])
            
        end
        function plot_cavity_radius(t, cavity_radius, param_vec, varargin)
            
            % setting optional parameter values, including parsing for
            % validation. parameters are: 
            % x_lim_end - setting the x_lim end value
            % x_units - setting the units for the x-axis
            default_x_lim_end = t(end);
            default_x_units = 'mins';
            expected_x_units = {'secs','mins','hrs','days'};

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            addParameter(p, 'x_units', default_x_units, ...
                @(x) any(validatestring(x,expected_x_units)));
            
            parse(p,varargin{:});
            
            switch p.Results.x_units
                case 'mins'
                    divisor = 60;
                case 'hrs'
                    divisor = 60*60;
                case 'days'
                    divisor = 24*60*60;
                otherwise
                    divisor = 1;
            end
            
            %plotting
            arrayLegend=cell(length(param_vec),1);
            for i = 1:length(param_vec)
                arrayLegend{i} = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
            end
            
            f = figure();
            ax = axes('Parent', f);

            cm = colormap('copper');
            
            [numColors,~] = size(cm);
            cm_divisor = floor(numColors/length(param_vec));
            
            [~, col] = size(cavity_radius);
            
            for i = 1:col
                plot(ax,t(:,i)/divisor, cavity_radius(:,i), 'LineWidth', 4, 'color', cm(cm_divisor*i,:));
                hold on
            end
            x_axis_label = sprintf('Time (%s)', p.Results.x_units);
            xlabel(ax,x_axis_label, 'FontSize', 32)
            ylabel(ax,'Cavity Radius (cm)','FontSize', 32)
            set(ax,'FontSize',28)
            legend(ax,arrayLegend)
            xlim(ax,[0, p.Results.x_lim_end/divisor])
            
        end
        function plot_concentration_profiles_cell(r, t, tvec, t_inds, C, param_vec, varargin)
            
            default_group_by = 'parameter';
            default_x_lim_end = r(end);
            default_ind_r_cavity = {};
            expected_group_by = {'parameter','time'};

            p = inputParser;

            addOptional(p, 'x_lim_end', default_x_lim_end);
            addOptional(p, 'ind_r_cavity', default_ind_r_cavity);
            addParameter(p, 'group_by', default_group_by, ...
                @(x) any(validatestring(x,expected_group_by)));
            
            parse(p,varargin{:});
            
            if strcmp(p.Results.group_by,'parameter') == 1
                
                for i = 1:length(param_vec)
                
                    titleLabel = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
                    
                    f = figure();
                    ax = axes('Parent', f);

                    cm = colormap('copper');
                    
                    [numColors,~] = size(cm);
                    cm_divisor = floor(numColors/length(tvec));
                    
                    conc = C{i};
                    
                    tvec_phi = t_inds{i};

                    for jj=1:length(tvec_phi)
                        plot(ax, r,conc(tvec_phi(jj),:),'LineWidth', 4, 'color', cm(cm_divisor*jj,:))
                        hold all
                    end

                    title(ax, titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                    ylabel(ax, 'Ethanol Concentration (g/ml)','FontSize',32,'FontWeight','Bold')
                    xlabel(ax, 'Radius (cm)','FontSize',32,'FontWeight','Bold')
                    set(ax,'FontSize',28)
                    
                    legendLabels = cell(length(t_inds),1);
                    for label = 1:length(legendLabels)
                    legendLabels{label} = sprintf('%d mins', t(t_inds(label),i)/60);
                    end
                    
                    legend(legendLabels)

                    xlim(ax, [0, p.Results.x_lim_end])
                    ylim(ax, [0, 1])

                end
                
            elseif strcmp(p.Results.group_by,'time') == 1
                
                for i = 1:length(tvec)
                    
                    titleLabel = sprintf('%d mins',tvec(i)/60);
                    
                    f = figure();
                    ax = axes('Parent', f);
                    
                    cm = colormap('copper');
                    
                    [numColors,~] = size(cm);
                    cm_divisor = floor(numColors/length(param_vec));
                    

                    for jj=1:length(param_vec)
                        conc = C{jj};
                        t_inds_vec = t_inds{jj};
                        plot(ax,r,conc(t_inds_vec(i),:),'LineWidth', 4, 'color', cm(cm_divisor*jj,:))
                        hold all
                        if numel(p.Results.ind_r_cavity) > 0
                            ind_r_cavity_phi = p.Results.ind_r_cavity{jj};
                            ind_r_cavity_at_t = ind_r_cavity_phi(t_inds_vec(i));
                            plot(ax, r(ind_r_cavity_at_t)*ones(100,1), linspace(0,1, 100), 'LineWidth',2, 'LineStyle', '--', 'color',cm(cm_divisor*jj,:))
                        end
                    end

                    title(ax, titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                    ylabel(ax, 'Ethanol Concentration (g/ml)','FontSize',32,'FontWeight','Bold')
                    xlabel(ax, 'Radius (cm)','FontSize',32,'FontWeight','Bold')
                    set(ax,'FontSize',28)
                    
                    if numel(p.Results.ind_r_cavity) ==  0
                        arrayLegend=cell(length(param_vec),1);
                        for param = 1:length(param_vec)
                            arrayLegend{param} = sprintf('%.0f%%w/w EC',(1-param_vec(param))*100);
                        end
                    else
                        arrayLegend={};
                        for param = 1:length(param_vec)
                            arrayLegend{end+1} =sprintf('%.0f%%w/w EC',(1-param_vec(param))*100);                
                            if param < length(param_vec)
                                arrayLegend{end+1} = '';
                            end
                        end
                        arrayLegend{end+1} = 'Rc';
                    end
                    legend(ax, arrayLegend)
              
                    xlim(ax, [0, p.Results.x_lim_end])
                    ylim(ax, [0, 1])

                end
                
            end
            
        end
        function plot_concentration_profiles(r, t, tvec, C, param_vec, varargin)
            
            default_group_by = 'parameter';
            default_x_lim_end = r(end);
            default_ind_r_cavity = [];
            expected_group_by = {'parameter','time'};

            p = inputParser;

            addOptional(p, 'x_lim_end', default_x_lim_end);
            addOptional(p, 'ind_r_cavity', default_ind_r_cavity);
            addParameter(p, 'group_by', default_group_by, ...
                @(x) any(validatestring(x,expected_group_by)));
            
            parse(p,varargin{:});
            
            if strcmp(p.Results.group_by,'parameter') == 1
                
                for i = 1:length(param_vec)
                
                    titleLabel = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
                    
                    figure(1+i)

                    cm = colormap('copper');
                    
                    [numColors,~] = size(cm);
                    cm_divisor = floor(numColors/length(tvec));
                    
                    conc = C{i};

                    for jj=1:length(tvec)
                        plot(r,conc(tvec(jj),:),'LineWidth', 4, 'color', cm(cm_divisor*jj,:))
                        hold all
                    end

                    title(titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                    ylabel('Ethanol Concentration (g/ml)','FontSize',32,'FontWeight','Bold')
                    xlabel('Radius (cm)','FontSize',32,'FontWeight','Bold')
                    set(gca,'FontSize',28)
                    
                    legendLabels = cell(length(tvec),1);
                    for label = 1:length(legendLabels)
                    legendLabels{label} = sprintf('%d mins', t(tvec(label),i)/60);
                    end
                    
                    legend(legendLabels)

                    xlim([0, p.Results.x_lim_end])
                    ylim([0, 1])

                end
                
            elseif strcmp(p.Results.group_by,'time') == 1
                
                for i = 1:length(tvec)
                    
                    titleLabel = sprintf('%d mins',t(tvec(i),1)/60);
                    
                    figure(1+i)

                    cm = colormap('copper');
                    
                    [numColors,~] = size(cm);
                    cm_divisor = floor(numColors/length(param_vec));
                    

                    for jj=1:length(param_vec)
                        conc = C{jj};
                        plot(r,conc(tvec(i),:),'LineWidth', 4, 'color', cm(cm_divisor*jj,:))
                        hold all
                        if numel(p.Results.ind_r_cavity) > 0
                            ind_r_cavity_at_t = p.Results.ind_r_cavity(tvec(i), jj);
                            plot(r(ind_r_cavity_at_t)*ones(100,1), linspace(0,1, 100), 'LineWidth',2, 'LineStyle', '--', 'color',cm(cm_divisor*jj,:))
                        end
                    end

                    title(titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                    ylabel('Ethanol Concentration (g/ml)','FontSize',32,'FontWeight','Bold')
                    xlabel('Radius (cm)','FontSize',32,'FontWeight','Bold')
                    set(gca,'FontSize',28)
                    
                    if numel(p.Results.ind_r_cavity) ==  0
                        arrayLegend=cell(length(param_vec),1);
                        for param = 1:length(param_vec)
                            arrayLegend{param} = sprintf('%.0f%%w/w EC',(1-param_vec(param))*100);
                        end
                    else
                        arrayLegend={};
                        for param = 1:length(param_vec)
                            arrayLegend{end+1} =sprintf('%.0f%%w/w EC',(1-param_vec(param))*100);                
                            if param < length(param_vec)
                                arrayLegend{end+1} = '';
                            end
                        end
                        arrayLegend{end+1} = 'Rc';
                    end
                    legend(arrayLegend)
              
                    xlim([0, p.Results.x_lim_end])
                    ylim([0, 1])

                end
                
            end
            
        end
        function plot_velocity_cell(t, tvec, r, ind_r0, ind_r_cavity, v, param_vec, varargin)
            
            default_x_lim_end = r(end);

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            
            parse(p,varargin{:});
            
            %figure()
            cm = colormap('copper');
            [numColors,~] = size(cm);
            
            cm_divisor = floor(numColors/length(param_vec));
            
            % v = cell -> 1xlength(phiVec)
            % v{i} = array -> length(r) x length(tvec)
            % want to plot: 1:length(tvec) plots of v vs. r for different
            % values of phi
            
            v_plot = cell(length(tvec),1);
            for i = 1:length(tvec)
                new_v = zeros(length(r), length(v));
                for j = 1:length(v)
                    v_temp = v{j};
                    new_v(:,j) = v_temp(:,i);
                end
                v_plot{i} = new_v;
            end
            
            for i = 1:length(tvec)
                figure()
                v_to_plot = v_plot{i};
                for j = 1:length(v)
                    plot(r, v_to_plot(:,j), 'color', cm(j*cm_divisor,:), 'LineWidth', 4)
                    hold all
                end
                xlabel('Radius (cm)', 'FontSize', 32)
                ylabel('Velocity (cm/s)','FontSize', 32)
                set(gca,'FontSize',28)
                xlim([0, p.Results.x_lim_end])
                title_text = sprintf('%d mins', t(tvec(i))/60);
                title(title_text, 'FontSize', 40)
                
                hold on
                plot(r(ind_r0)*ones(100,1), linspace(min(min(v_to_plot)),max(max(v_to_plot)),100),'LineStyle', '--', 'LineWidth', 2, 'color', [0.75, 0.75, 0.75])
                %plot(r(ind_r_cavity(1))*ones(100,1), linspace(min(min(v_to_plot)),max(max(v_to_plot)),100),'LineStyle', '--', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4])
                
                legendLabels = cell(length(param_vec)+1,1);
                for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                    legendLabels{label} = sprintf('%.0f%%w/w EC',(1-param_vec(label))*100);
                end
                legendLabels{end} = 'r_{0}';     
                %legendLabels{end} = 'r_{c}';             

                legend(legendLabels)

            end
            
        end
        function plot_velocity(t, tvec, r, ind_r0, ind_r_cavity, v, param_vec, varargin)
            
            default_x_lim_end = r(end);

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            
            parse(p,varargin{:});
            
            %figure()
            cm = colormap('copper');
            [numColors,~] = size(cm);
            
            cm_divisor = floor(numColors/length(param_vec));
            
            % v = cell -> 1xlength(phiVec)
            % v{i} = array -> length(r) x length(tvec)
            % want to plot: 1:length(tvec) plots of v vs. r for different
            % values of phi
            
            v_plot = cell(length(tvec),1);
            for i = 1:length(tvec)
                new_v = zeros(length(r), length(v));
                for j = 1:length(v)
                    v_temp = v{j};
                    new_v(:,j) = v_temp(:,i);
                end
                v_plot{i} = new_v;
            end
            
            for i = 1:length(tvec)
                figure()
                v_to_plot = v_plot{i};
                for j = 1:length(v)
                    plot(r, v_to_plot(:,j), 'color', cm(j*cm_divisor,:), 'LineWidth', 4)
                    hold all
                end
                xlabel('Radius (cm)', 'FontSize', 32)
                ylabel('Velocity (cm/s)','FontSize', 32)
                set(gca,'FontSize',28)
                xlim([0, p.Results.x_lim_end])
                title_text = sprintf('%d mins', t(tvec(i))/60);
                title(title_text, 'FontSize', 40)
                
                hold on
                plot(r(ind_r0)*ones(100,1), linspace(min(min(v_to_plot)),max(max(v_to_plot)),100),'LineStyle', '--', 'LineWidth', 2, 'color', [0.75, 0.75, 0.75])
                %plot(r(ind_r_cavity(1))*ones(100,1), linspace(min(min(v_to_plot)),max(max(v_to_plot)),100),'LineStyle', '--', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4])
                
                legendLabels = cell(length(param_vec)+1,1);
                for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                    legendLabels{label} = sprintf('%.0f%%w/w EC',(1-param_vec(label))*100);
                end
                legendLabels{end} = 'r_{0}';     
                %legendLabels{end} = 'r_{c}';             

                legend(legendLabels)

            end
            
        end
        function plot_mass_conservation_cell(t, a_vec, r, ind_r0, ind_r_cavity, C, param_vec, time_pts, tvec_plot, varargin)
            
            plot_line_plots = 0;

            temp = t{1};
            default_t_units = 'mins';
            default_x_lim_end = temp(end);
            expected_t_units = {'mins','hrs','days'};

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            addParameter(p, 't_units', default_t_units, ...
                @(x) any(validatestring(x,expected_t_units)));
            
            parse(p,varargin{:});
            
            time_divisor = 60;
            time_label = 'mins';
            
            if strcmp(p.Results.t_units,'hrs') == 1
                time_divisor = 3600;
                time_label = 'hrs';
            elseif strcmp(p.Results.t_units,'days') == 1
                time_divisor = 24*3600;
                time_label = 'days';
            end
            
            %figure()
            cm = colormap('copper');
            [numColors,~] = size(cm);
            
            cm_divisor = floor(numColors/length(param_vec));
            
            %mass_inside_cavity = zeros(length(t), length(param_vec));
            %mass_inside_tumor = zeros(length(t), length(param_vec));
            %mass_inside_tissue = zeros(length(t), length(param_vec));
            %mass_inside_cavity = cell(length(param_vec),1);
            %mass_inside_tumor = cell(length(param_vec),1);
            %mass_inside_tissue = cell(length(param_vec),1);

            
            legendLabels = cell(length(param_vec),1);
            
            for i = 1:length(param_vec)
                t_vec = t{i};
                tvec_plot_vec = tvec_plot{i};
                mass_inside_cavity_vec = zeros(length(t_vec),1);
                mass_inside_tumor_vec = zeros(length(t_vec),1);
                mass_inside_tissue_vec = zeros(length(t_vec),1);
                total_mass = zeros(length(t_vec),1);

                
                ind_r_cavity_vec = ind_r_cavity{i};
                a_phi = a_vec{i};
                conc = C{i};

                for j = 1:length(t_vec)
                    bound_cavity = ind_r_cavity_vec(j);    
                                              %4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1)
                    mass_inside_cavity_vec(j) = 4 * pi * trapz(r(1:bound_cavity), (r(1:bound_cavity).^2).*conc(j,1:bound_cavity), 2);
                    mass_inside_tumor_vec(j) = 4 * pi * trapz(r(bound_cavity+1:a_phi(j)), (r(bound_cavity+1:a_phi(j)).^2).*conc(j,bound_cavity+1:a_phi(j)), 2);
                    mass_inside_tissue_vec(j) = 4 * pi * trapz(r(a_phi(j)+1:end), (r(a_phi(j)+1:end).^2).*conc(j,a_phi(j)+1:end), 2);
                    total_mass(j) = mass_inside_cavity_vec(j) + mass_inside_tumor_vec(j) + mass_inside_tissue_vec(j);

                end

                stack = zeros(length(time_pts),3);
                xlabel_values = cell(length(time_pts),1);

                for k=1:length(time_pts)
                    stack(k,1) = mass_inside_cavity_vec(tvec_plot_vec(k))/total_mass(tvec_plot_vec(k));
                    stack(k,2) = mass_inside_tumor_vec(tvec_plot_vec(k))/total_mass(tvec_plot_vec(k));
                    stack(k,3) = -1*mass_inside_tissue_vec(tvec_plot_vec(k))/total_mass(tvec_plot_vec(k));
                    xlabel_values{k} = sprintf('%.0f',time_pts(k)/time_divisor);
                end
                
                f = figure;
                ax = axes('Parent',f);
                bar(ax, stack, 'stacked')
                ylabel(ax, 'Normalized Ethanol Fraction','FontSize',12,'FontWeight','Bold')
                xlabel(ax, sprintf('Time (%s)', time_label),'FontSize',12,'FontWeight','Bold')
                title(ax, sprintf('%.0f%%w/w EC',(1-param_vec(i))*100),'FontSize',14,'FontWeight','Bold')
                xticklabels(ax, xlabel_values)
                set(ax,'FontSize', 28)
                ylim([-0.4 1])
                legend(ax, {'Cavity','Tumour','Tissue'})

                
                if plot_line_plots == 1
                    if i == 1
                        f1 = figure;
                        f2 = figure;
                        f3 = figure;
                        
                        ax1 = axes('Parent', f1);
                        ax2 = axes('Parent', f2);
                        ax3 = axes('Parent', f3);
                        
                        hold(ax1,'on');
                        hold(ax2,'on');
                        hold(ax3,'on');
                        
                        set(ax1,'FontSize',28);
                        set(ax2,'FontSize',28);
                        set(ax3,'FontSize',28);
    
                    end
                    plot(ax1, t_vec/time_divisor, mass_inside_cavity_vec, 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                    plot(ax2, t_vec/time_divisor, mass_inside_tumor_vec, 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                    plot(ax3, t_vec/time_divisor, mass_inside_tissue_vec, 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                    
                    legendLabels{i} = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
                end
            end
            
            if plot_line_plots == 1
                title(ax1, 'Mass OH inside cavity', 'FontSize', 40, 'FontWeight', 'Bold')
                title(ax2, 'Mass OH inside tumor', 'FontSize', 40, 'FontWeight', 'Bold')
                title(ax3, 'Mass OH inside tissue', 'FontSize', 40, 'FontWeight', 'Bold')
                xlabel_text = sprintf('Time (%s)', time_label);
                xlabel([ax1,ax2,ax3], xlabel_text, 'FontSize', 32)
                ylabel([ax1,ax2,ax3], 'Mass Ethanol', 'FontSize', 32)
                legend(ax1, legendLabels)
                legend(ax2, legendLabels)
                legend(ax3, legendLabels)
                xlim(ax1, [0, p.Results.x_lim_end/time_divisor])
                xlim(ax2, [0, p.Results.x_lim_end/time_divisor])
                xlim(ax3, [0, p.Results.x_lim_end/time_divisor])
            end
                                           
        end
        function plot_mass_conservation(t, a_vec, r, ind_r0, ind_r_cavity, C, param_vec, varargin)
            
            default_t_units = 'mins';
            default_x_lim_end = t(end);
            expected_t_units = {'mins','hrs','days'};

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            addParameter(p, 't_units', default_t_units, ...
                @(x) any(validatestring(x,expected_t_units)));
            
            parse(p,varargin{:});
            
            time_divisor = 60;
            time_label = 'mins';
            
            if strcmp(p.Results.t_units,'hrs') == 1
                time_divisor = 3600;
                time_label = 'hrs';
            elseif strcmp(p.Results.t_units,'days') == 1
                time_divisor = 24*3600;
                time_label = 'days';
            end
            
            %figure()
            cm = colormap('copper');
            [numColors,~] = size(cm);
            
            cm_divisor = floor(numColors/length(param_vec));
            
            mass_inside_cavity = zeros(length(t), length(param_vec));
            mass_inside_tumor = zeros(length(t), length(param_vec));
            mass_inside_tissue = zeros(length(t), length(param_vec));
            
            legendLabels = cell(length(param_vec),1);
            
            for i = 1:length(param_vec)
                for j = 1:length(t)
                    bound_cavity = ind_r_cavity(j,i);
                    conc = C{i};
                                              %4 * pi * trapz(r(1:ind_Y_r_cavity), (r(1:ind_Y_r_cavity)'.^2) .* Y(1:ind_Y_r_cavity,:), 1)
                    mass_inside_cavity(j,i) = 4 * pi * trapz(r(1:bound_cavity), (r(1:bound_cavity).^2).*conc(j,1:bound_cavity), 2);
                    mass_inside_tumor(j,i) = 4 * pi * trapz(r(bound_cavity:a_vec(j,i)), (r(bound_cavity:a_vec(j,i)).^2).*conc(j,bound_cavity:a_vec(j,i)), 2);
                    mass_inside_tissue(j,i) = 4 * pi * trapz(r(a_vec(j,i):end), (r(a_vec(j,i):end).^2).*conc(j,a_vec(j,i):end), 2);
                end
                
                if i == 1
                    f1 = figure;
                    f2 = figure;
                    f3 = figure;
                    
                    ax1 = axes('Parent', f1);
                    ax2 = axes('Parent', f2);
                    ax3 = axes('Parent', f3);
                    
                    hold(ax1,'on');
                    hold(ax2,'on');
                    hold(ax3,'on');
                    
                    set(ax1,'FontSize',28);
                    set(ax2,'FontSize',28);
                    set(ax3,'FontSize',28);

                end
                plot(ax1, t/time_divisor, mass_inside_cavity(:,i), 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                plot(ax2, t/time_divisor, mass_inside_tumor(:,i), 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                plot(ax3, t/time_divisor, mass_inside_tissue(:,i), 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                
                legendLabels{i} = sprintf('%.0f%%w/w EC',(1-param_vec(i))*100);
            end
            
            title(ax1, 'Mass OH inside cavity', 'FontSize', 40, 'FontWeight', 'Bold')
            title(ax2, 'Mass OH inside tumor', 'FontSize', 40, 'FontWeight', 'Bold')
            title(ax3, 'Mass OH inside tissue', 'FontSize', 40, 'FontWeight', 'Bold')
            xlabel_text = sprintf('Time (%s)', time_label);
            xlabel([ax1,ax2,ax3], xlabel_text, 'FontSize', 32)
            ylabel([ax1,ax2,ax3], 'Mass Ethanol', 'FontSize', 32)
            legend(ax1, legendLabels)
            legend(ax2, legendLabels)
            legend(ax3, legendLabels)
            xlim(ax1, [0, p.Results.x_lim_end/time_divisor])
            xlim(ax2, [0, p.Results.x_lim_end/time_divisor])
            xlim(ax3, [0, p.Results.x_lim_end/time_divisor])
            
                                    
%             mass_injected_tot = (C0_OH*Vol).*ones(1,length(t_OH));
% 
%             bound = find(r >= floor(a(1)), 1);
% 
%             mass_inside_oh = 4*pi*trapz((r(1:bound)'.^2).*C_OH(1:bound,:),1)/(length(r(1:bound))-1);
%             mass_outside_oh = 4*pi*trapz((r(bound+1:end)'.^2).*C_OH(bound+1:end,:),1)/(length(r(1:bound))-1);
%             mass_injected_oh = mass_inside_oh+mass_outside_oh;

%             figure()
%             plot(t_OH/60,mass_inside_oh,'LineWidth',2)
%             hold on
%             plot(t_OH/60,mass_outside_oh,'r','LineWidth',2)
%             plot(t_OH/60,mass_injected_oh,'k--','LineWidth',2)
%             plot(t_OH/60,mass_injected_tot,'k','LineWidth',2)
% 
%             ylabel('Ethanol Mass (mg)','FontSize',18,'FontWeight','Bold')
%             xlabel('Time (mins)','FontSize',18,'FontWeight','Bold')
%             title('Total Mass of Ethanol Inside and Outside Tumor','FontSize',18,'FontWeight','Bold')
%             XT = get(gca,'XTick');
%             set(gca,'FontSize',16)
%             YT = get(gca,'YTick');
%             set(gca,'FontSize',16)
%             legend('Inside tumor', 'Outside Tumor', 'Inside + outside', 'Total injected')
        end
        function plot_ethanol_cloud_cell(t, t_vec, t_inds, C, r, c_threshold, param_vec, varargin)
            ethanol_cloud = zeros(length(t_vec), length(param_vec));
            %EC_conc = zeros(length(param_vec), 1);
            
            xlabels = cell(length(param_vec),1);
            for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                xlabels{label} = sprintf('%.0f%%',(1-param_vec(label))*100);
            end
            
            X = categorical(xlabels);
            X = reordercats(X, xlabels);
            
            for i = 1:length(t_vec)
                for j = 1:length(param_vec)
                   t_inds_param = t_inds{j};
                   conc = C{j};
                   r_threshold_indx = find(conc(t_inds_param(i),:) <= c_threshold, 1);
                   r_cloud = r(r_threshold_indx);
                   ethanol_cloud(i,j) = (4/3)*pi*r_cloud.^3;
                end
               %EC_conc(i) = (1-(param_vec(i)))*100;
                f = figure;
                ax = axes('Parent',f);
                bar(ax, X, ethanol_cloud*1000)
                set(ax,'FontSize', 28)
                ylabel(ax, 'Ethanol cloud volume (uL)', 'FontSize', 32)
                xlabel(ax, 'EC concentration (%w/w)', 'FontSize',32)
            end
            

         end
         function plot_ethanol_cloud(t, t_vec, C, r, c_threshold, param_vec, varargin)
            ethanol_cloud = zeros(length(t_vec), length(param_vec));
            %EC_conc = zeros(length(param_vec), 1);
            
            xlabels = cell(length(param_vec),1);
            for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                xlabels{label} = sprintf('%.0f%%',(1-param_vec(label))*100);
            end
            
            X = categorical(xlabels);
            X = reordercats(X, xlabels);
            
            for i = 1:length(t_vec)
                for j = 1:length(param_vec)
                   conc = C{j};
                   r_threshold_indx = find(conc(t_vec(i),:) <= c_threshold, 1);
                   r_cloud = r(r_threshold_indx);
                   ethanol_cloud(i,j) = (4/3)*pi*r_cloud.^3;
                end
               %EC_conc(i) = (1-(param_vec(i)))*100;
                f = figure;
                ax = axes('Parent',f);
                bar(ax, X, ethanol_cloud*1000)
                set(ax,'FontSize', 28)
                ylabel(ax, 'Ethanol cloud volume (uL)', 'FontSize', 32)
                xlabel(ax, 'EC concentration (%w/w)', 'FontSize',32)
            end
            

         end
         function plot_cavity_volume_cell(t, a, r, cavity_radius, param_vec, varargin)
            %cavity_volume = zeros(length(t), length(param_vec));
            %tumor_volume = zeros(length(t), length(param_vec));
            
            cavity_volume = cell(length(param_vec),1);
            tumor_volume = cell(length(param_vec),1);

            
            xlabels = cell(length(param_vec),1);
            for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                xlabels{label} = sprintf('%.0f%%',(1-param_vec(label))*100);
            end
            
            X = categorical(xlabels);
            X = reordercats(X, xlabels);
            t_pt = 5*60;
            t_inds = zeros(length(param_vec),1);
            end_cavity_volumes = zeros(length(param_vec),1);
            
            for i = 1:length(param_vec)
                cavity_radius_vec = cavity_radius{i};
                t_vec = t{i};
                a_vec = a{i};
                cavity_volume_vec = (4/3)*pi*cavity_radius_vec.^3;
                cavity_volume{i} = cavity_volume_vec;
                tumor_volume{i} = (4/3)*pi*a_vec.^3-cavity_volume_vec;
                t_inds(i) = find(t_vec >= t_pt, 1);
                end_cavity_volumes(i) = cavity_volume_vec(t_inds(i),:);

                %EC_conc(i) = (1-(param_vec(i)))*100;
            end
            
            %t_ind = find(t>=5*60,1);
            
            f = figure;
            ax = axes('Parent',f);
            bar(ax, X, end_cavity_volumes*1000)
            set(ax,'FontSize', 28)
            ylabel(ax, 'Cavity volume (uL)', 'FontSize', 32)
            xlabel(ax, 'EC concentration (%w/w)', 'FontSize',32)
            
%             xlabel('Time (mins)', 'FontSize', 32)
            
%             f = figure;
%             ax1 = axes('Parent', f);
%             cm = colormap('copper');
%             [numColors,~] = size(cm);
%             cm_divisor = floor(numColors/length(param_vec));
%             
%             legendLabels = cell(length(phi_val)*2,1);
%             
%             for i = 1:length(param_vec)
%                 plot(ax1, tt/60, a(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%                 hold on
%                 plot(ax1, tt/60, cavity_radius(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
%                 legendLabels{i} = sprintf('%.0f%%w/w EC, tumor',(1-phiVec(i))*100);
%             end
%             for label = length(param_vec)+1:2*length(param_vec)
%                 legendLabels{label} = sprintf('%.0f%%w/w EC, cavity',(1-param_vec(floor(label/2)))*100);
%             end
%             
%             legend(ax1, legendLabels)
%             set(ax1,'FontSize', 28)
%             xlim([0, 6])
%             ylabel('Radius (cm)', 'FontSize', 32)
%             xlabel('Time (mins)', 'FontSize', 32)
%             
%             f = figure;
%             ax2 = axes('Parent',f);
%             for i = 1:length(param_vec)
%                 plot(ax2, tt/60, cavity_volume(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%                 hold on
%                 plot(ax2, tt/60, tumor_volume(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
%             end
%             
%             legend(ax2, legendLabels)
%             set(ax2,'FontSize', 28)
%             xlim([0, 6])
%             ylabel('Volume (cm^{3})', 'FontSize', 32)
%             xlabel('Time (mins)', 'FontSize', 32)

        end
        function plot_cavity_volume(t, a, r, cavity_radius, param_vec, varargin)
            cavity_volume = zeros(length(t), length(param_vec));
            tumor_volume = zeros(length(t), length(param_vec));
            %EC_conc = zeros(length(param_vec), 1);
            
            xlabels = cell(length(param_vec),1);
            for label = 1:length(param_vec)
                    %legendLabels{label} = sprintf('%d mins', t(tvec(label))/60);
                xlabels{label} = sprintf('%.0f%%',(1-param_vec(label))*100);
            end
            
            X = categorical(xlabels);
            X = reordercats(X, xlabels);
            
            for i = 1:length(param_vec)
               cavity_volume(:,i) = (4/3)*pi*cavity_radius(:,i).^3;
               tumor_volume(:,i) = (4/3)*pi*a(:,i).^3-cavity_volume(:,i);
               %EC_conc(i) = (1-(param_vec(i)))*100;
            end
            
            t_ind = find(t>=5*60,1);
            end_cavity_volume = cavity_volume(t_ind,:);
            
            f = figure;
            ax = axes('Parent',f);
            bar(ax, X, end_cavity_volume*1000)
            set(ax,'FontSize', 28)
            ylabel(ax, 'Cavity volume (uL)', 'FontSize', 32)
            xlabel(ax, 'EC concentration (%w/w)', 'FontSize',32)
            
%             xlabel('Time (mins)', 'FontSize', 32)
            
%             f = figure;
%             ax1 = axes('Parent', f);
%             cm = colormap('copper');
%             [numColors,~] = size(cm);
%             cm_divisor = floor(numColors/length(param_vec));
%             
%             legendLabels = cell(length(phi_val)*2,1);
%             
%             for i = 1:length(param_vec)
%                 plot(ax1, tt/60, a(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%                 hold on
%                 plot(ax1, tt/60, cavity_radius(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
%                 legendLabels{i} = sprintf('%.0f%%w/w EC, tumor',(1-phiVec(i))*100);
%             end
%             for label = length(param_vec)+1:2*length(param_vec)
%                 legendLabels{label} = sprintf('%.0f%%w/w EC, cavity',(1-param_vec(floor(label/2)))*100);
%             end
%             
%             legend(ax1, legendLabels)
%             set(ax1,'FontSize', 28)
%             xlim([0, 6])
%             ylabel('Radius (cm)', 'FontSize', 32)
%             xlabel('Time (mins)', 'FontSize', 32)
%             
%             f = figure;
%             ax2 = axes('Parent',f);
%             for i = 1:length(param_vec)
%                 plot(ax2, tt/60, cavity_volume(:,i), 'LineWidth', 4, 'color', cm(i*cm_divisor,:))
%                 hold on
%                 plot(ax2, tt/60, tumor_volume(:,i), 'LineWidth', 4, 'LineStyle', '--', 'color', cm(i*cm_divisor,:))
%             end
%             
%             legend(ax2, legendLabels)
%             set(ax2,'FontSize', 28)
%             xlim([0, 6])
%             ylabel('Volume (cm^{3})', 'FontSize', 32)
%             xlabel('Time (mins)', 'FontSize', 32)

        end
        function plot_pc(t, pc, pcrit, param_vec, Vol, Q_needle, varargin)
            
            default_t_units = 'mins';
            default_x_lim_end = t(end);
            expected_t_units = {'mins','hrs','days'};

            p = inputParser;
            addOptional(p, 'x_lim_end', default_x_lim_end);
            addParameter(p, 't_units', default_t_units, ...
                @(x) any(validatestring(x,expected_t_units)));
            
            parse(p,varargin{:});
            
            time_divisor = 60;
            time_label = 'mins';
            
            if strcmp(p.Results.t_units,'hrs') == 1
                time_divisor = 3600;
                time_label = 'hrs';
            elseif strcmp(p.Results.t_units,'days') == 1
                time_divisor = 24*3600;
                time_label = 'days';
            end
            
            legendLabels = cell(length(param_vec)+1,1);
            legendLabels{1} = sprintf('Pcrit');
            for label = 2:length(param_vec)+1
                legendLabels{label} = sprintf('%.0f%% ECE',(1-param_vec(label-1))*100);
            end
            
            f = figure;
            ax = axes('Parent',f);
            
            cm = colormap('copper');
            [numColors,~] = size(cm);
            
            cm_divisor = floor(numColors/length(param_vec));
            
            tt_plot = t{1};
            pcrit_plot = pcrit{1};
            plot(ax, tt_plot, pcrit_plot, 'color', [0,0,0], 'LineWidth', 2, 'LineStyle', '--')
            hold on

            for i = 1:length(param_vec)
                pc_plot = pc{i};
                tt_plot = t{i};
                plot(ax, tt_plot/time_divisor, pc_plot, 'color', cm(i*cm_divisor,:), 'LineWidth', 4)
                hold on
            end
            
            xlabel_text = sprintf('Time (%s)', time_label);
            set(ax,'FontSize', 28)
            ylabel(ax, 'Pressure (Pa)', 'FontSize', 32)
            xlabel(ax, xlabel_text, 'FontSize',32)
            xlim([0, floor((Vol/Q_needle)/(time_divisor))])
            legend(ax, legendLabels)

            
        end
        function plot_strain()
            
        end
    end
end