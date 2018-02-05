%% Convert SubFREHD-C output files into .mat format
% Zhi Li 20170511

% path of the output data files
filepath = '/Users/lichee/Documents/Research/SubFREHD-C/WeirV5noinflow/output_10x10sub/';

% saved .mat file name
savefile = 1;
savename = 'SubgridTest_WeirV5noinflow_10x10sub.mat';

% dimension and time settings of the computation domain
% setting.Nx = 331;
% setting.Ny = 482;
% setting.Nx = 185;
% setting.Ny = 235;
% setting.Nx = 150;
% setting.Ny = 200;
% setting.Nx = 600;
% setting.Ny = 800;
setting.Nx = 60;
setting.Ny = 80;
% setting.Nx = 560;
% setting.Ny = 1408;
% setting.Nx = 280;
% setting.Ny = 704;
% setting.Nx = 140;
% setting.Ny = 352;
% setting.Nx = 70;
% setting.Ny = 176;
% setting.Nx = 35;
% setting.Ny = 88;
% setting.Nx = 333;
% setting.Ny = 485;

setting.dt = 1;
setting.t0 = '2013-06-10';
setting.subgrid = 1;
setting.saveCD = 1;

% index of file numbers
% nuecesUpChannelBay dt=1
% fileNumber = [0:1800:259200];
% fileNumber = [0:1800:72000];
% nuecesUpChannelBay dt=2
% fileNumber = [0:900:64800];
% nuecesUpChannelBay dt=5
% fileNumber = [0:360:25920];
% nuecesUpChannelBay dt=10
% fileNumber = [0:180:12960];
% Weir and Stair
fileNumber = [0:360:28800];
% fileNumber = [43200:360:86400];
% NDHM
% fileNumber = [0:720:15840];
Lt = length(fileNumber);

% initialize data arrays
data.uu = zeros(setting.Nx, setting.Ny, Lt);
data.vv = zeros(setting.Nx, setting.Ny, Lt);
data.surface = zeros(setting.Nx, setting.Ny, Lt);
data.depth = zeros(setting.Nx, setting.Ny, Lt);
data.scalar = zeros(setting.Nx, setting.Ny, Lt);
data.index = fileNumber;
data.time = (fileNumber .* setting.dt ./ 86400) + datenum(setting.t0);
if setting.subgrid == 1
    data.N = zeros(setting.Nx, setting.Ny, Lt);
    data.O = zeros(setting.Nx, setting.Ny, Lt);
    data.V = zeros(setting.Nx, setting.Ny, Lt);
    data.Z = zeros(setting.Nx, setting.Ny, Lt);
end
if setting.saveCD == 1
    data.CDXP = zeros(setting.Nx, setting.Ny, Lt);
    data.CDYP = zeros(setting.Nx, setting.Ny, Lt);
end

% saving data files
for jj = 1:Lt
    ii = fileNumber(jj);
    % depth
    fid = fopen([filepath,'depth_',num2str(ii),'.dat'],'r');
    temp = fscanf(fid, '%f');
    fclose(fid);
    aa = temp <= 0;
    temp(aa) = NaN;
    data.depth(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    % u velocity
    fid = fopen([filepath,'uVelocity_',num2str(ii),'.dat'],'r');
    temp = fscanf(fid, '%f');
    fclose(fid);
    temp(aa) = NaN;
    data.uu(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    % v velocity
    fid = fopen([filepath,'vVelocity_',num2str(ii),'.dat'],'r');
    temp = fscanf(fid, '%f');
    fclose(fid);
    temp(aa) = NaN;
    data.vv(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    % surface elevation
    fid = fopen([filepath,'surfaceZ_',num2str(ii),'.dat'],'r');
    temp = fscanf(fid, '%f');
    fclose(fid);
    temp(aa) = NaN;
    data.surface(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    % scalar
    fid = fopen([filepath,'scalar_',num2str(ii),'.dat'],'r');
    temp = fscanf(fid, '%f');
    fclose(fid);
    temp(aa) = NaN;
    data.scalar(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    % subgrid variables
    if setting.subgrid == 1
        fid = fopen([filepath,'SubNx_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        data.N(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
        fid = fopen([filepath,'SubOy_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        data.O(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
        fid = fopen([filepath,'SubV_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        data.V(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
        fid = fopen([filepath,'SubZ_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        data.Z(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    end
    % bottom drag coefficient
    if setting.saveCD == 1
        fid = fopen([filepath,'cdxp_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        temp(aa) = NaN;
        data.CDXP(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
        fid = fopen([filepath,'cdyp_',num2str(ii),'.dat'],'r');
        temp = fscanf(fid, '%f');
        fclose(fid);
        temp(aa) = NaN;
        data.CDYP(:,:,jj) = reshape(temp, setting.Nx, setting.Ny);
    end
end

% save datafile
if savefile == 1
    save(savename, 'data', '-v7.3');
end