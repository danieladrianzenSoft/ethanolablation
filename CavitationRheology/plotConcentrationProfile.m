function plotConcentrationProfile(r, tvec1, C_OH)

    figure()

    cm = colormap('copper');

    for jj=1:length(tvec1)

    plot(r,C_OH(:,tvec1(jj)),'color',cm(jj*20,:),'LineWidth',2)
    hold all

    end

    ylabel('Ethanol Concentration (normalized)','FontSize',18,'FontWeight','Bold')
    xlabel('Radius (cm)','FontSize',18,'FontWeight','Bold')
    title('Normalized ethanol concentration vs. radius at different times','FontSize',18,'FontWeight','Bold')
    legend('0mins','1min','2mins','4mins','6mins','8mins','12mins','20mins','30mins')
    XT = get(gca,'XTick');
    set(gca,'FontSize',16)
    YT = get(gca,'YTick');
    set(gca,'FontSize',16)
    xlim([0,0.8])

end
