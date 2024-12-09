load(lgd_color_txt)

m_proj('robinson','lon',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
[long,lati] = meshgrid(lon,lat);
m_pcolor(long,lati,field');
shading flat;
colormap(color);
hold on
if(plt_proxy)
    h=m_scatter(proxy_lon,proxy_lat,10,proxy_SNR,'LineWidth',0.9,'MarkerEdgeColor','k');
end
m_coast();
% m_grid('xtick',[],'ytick',[],'BackgroundColor',[0.7 0.7 0.7]);
m_grid('xtick',6,'tickdir','out','BackgroundColor',[0.7 0.7 0.7]);
c=colorbar('southoutside','fontsize',12);
caxis([left_lgd,right_lgd]);
c.Label.String = label_txt;
c.Label.FontSize = 15;
% ylabel(title_txt,'FontSize',12)
title(title_txt)
hold off