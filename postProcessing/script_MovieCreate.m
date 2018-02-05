%% Create movie for the SubFREHD-C outputs
% Zhi Li 20170522

% name of the data file to read
fname = 'SubgridTest_WeirV5noinflow_10x10sub.mat';
% variable of the movie
dataname = 'scalar';

% other movie settings
setting.figposition = [100 100 800 600];%[100 100 550 600];
setting.timeLoc = [0.12 0.2 0.3 0.1];
setting.colorRange = [-2 27];
setting.cmap = colormap(jet);
setting.fs = 16;
setting.title = 'Salinity [psu]';

% load data 
load(fname);
Mdata = data.(dataname);
Nt = length(data.time);

% aux data
cmin = setting.colorRange(1);
cmax = setting.colorRange(2);
cdiff = cmax - cmin;
cmap = setting.cmap;
xVec = 1:size(Mdata,1);
yVec = size(Mdata,2):-1:1;
% create the movie object
fig = figure(1);
iframe = 0;
movObj = VideoWriter([dataname,'.avi']);
open(movObj);
for ii = 1:Nt
    iframe = iframe + 1;
    set(gcf,'Position',setting.figposition);
    piece = Mdata(:,:,ii);
    imagesc(xVec, yVec, piece);
    set(gca,'Ydir','normal');
    caxis([cmin - cdiff/64, cmax]);
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    title(setting.title,'FontSize',setting.fs);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    annotation(fig,'textbox',setting.timeLoc,'string',...
        {datestr(data.time(ii))},'FontSize',setting.fs,...
        'Color',[0 0 0],'EdgeColor',[1 1 1],...
        'FontWeight','bold','LineStyle','none');
    cb = colorbar('Location','eastoutside');
    set(cb,'FontSize',setting.fs);
    ylabel(cb,setting.title,'FontSize',setting.fs);
    % capture the frame
    writeVideo(movObj, getframe(gcf));
    if ii < Nt
        clf
    end
end
close(movObj);


