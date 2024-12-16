load(lgd_color_txt)
Graph{1,1}=rgb('Gold');        Graph{1,2}= 'h'; % bivalve
Graph{2,1}=rgb('DarkKhaki');   Graph{2,2}= 'h'; % borehole
Graph{3,1}=rgb('DarkOrange');  Graph{3,2}= 'o'; % coral
Graph{4,1}=rgb('Black');       Graph{4,2}= 'p'; % documents
Graph{5,1}=rgb('LightSkyBlue');Graph{5,2}= 'd'; % glacier ice
Graph{6,1}=rgb('DeepSkyBlue'); Graph{6,2}= '*'; % hybrid
Graph{7,1}=rgb('RoyalBlue');   Graph{7,2}= 's'; % lake sediment
Graph{8,1}=rgb('SaddleBrown'); Graph{8,2}= 's'; % marine sediment
Graph{9,1}=rgb('Red');         Graph{9,2}= 'o'; % sclerosponge
Graph{10,1}=rgb('DeepPink');   Graph{10,2}= 'd'; % speleothem
Graph{11,1}=rgb('LimeGreen');  Graph{11,2}= '^'; % tree

if(~isempty(find(cellfun(@(x) contains(x, 'coral'), Ptype{1}(filter_space)))))
    p_code = 3;
elseif(~isempty(find(cellfun(@(x) contains(x, 'tree'), Ptype{1}(filter_space)))))
    p_code = 11;
elseif(~isempty(find(cellfun(@(x) contains(x, 'ice'), Ptype{1}(filter_space)))))
    p_code = 5;
elseif(~isempty(find(cellfun(@(x) contains(x, 'lake'), Ptype{1}(filter_space)))))
    p_code = 7;
elseif(~isempty(find(cellfun(@(x) contains(x, 'bivalve'), Ptype{1}(filter_space)))))
    p_code = 1;
end

m_proj('robinson','lon',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
[long,lati] = meshgrid(lon,lat);
m_pcolor(long,lati,field');
shading flat;
colormap(color);
hold on
if(plt_proxy)
    m_line(proxy_lon, proxy_lat, 'marker',Graph{p_code,2},'MarkerEdgeColor','black','MarkerFaceColor',Graph{p_code,1},'linewidth',[2],'MarkerSize',[8],'linestyle','none');    
end
h = m_line(x_lon, x_lat, 'marker','o','MarkerEdgeColor','black','linewidth',[2],'MarkerSize',[6],'linestyle','none');
m_coast();
m_grid('xtick',6,'tickdir','out','BackgroundColor',[0.7 0.7 0.7]);
c = colorbar('southoutside','fontsize',12);
caxis([left_lgd,right_lgd]);
c.Label.String = label_txt;
c.Label.FontSize = 15;
title(title_txt)
hold off