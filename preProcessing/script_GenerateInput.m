%% This script generates input files for the SubFREHD-C
%% This script is customized for the Nueces Delta model
% Zhi Li 20170608

% path(path, 'C:\Users\leach\Documents\Research\NuecesDelta\');
% path(path, '~/Documents/Research/SubFREHD-C/input_Upper20x20NWD/');
fname = 'nuecesUCBV2_2x2.mat';
fnamePump = 'dailypump_const00.mat';
fnameic = 'scalarIC_nuecesUCBV22x2,mat';
fnameTide = 'tide_sinusoidal_T24R06.mat';
fnamesubA = 'subArea_1x1to10x10StairV3.mat';
fnamesubB = 'subBath_1x1to10x10StairV3.mat';

% files to be generated
ibathymetry = 0;
iedge = 0;
itide = 0;
iinflow = 1;
iwind = 0;
isalinity = 0;
iscalaric = 0;
isubgrid = 0;


% date number and the corresponding date string in C format
% tNum = 1325394000;% - 7200;
tNum = 1325397600;
tStr = '2012-01-01';
tNumM = datenum(tStr);

%% Generate bathymetry file
if ibathymetry == 1
    load(fname);
    data = bathCenter;
    dim = [size(data,1) size(data,2)];
    N = dim(1)*dim(2);
    aa = (isnan(data) | (data < -100));
    data(aa) = 10;
    bath1D = reshape(data, N, 1);
    fid = fopen('bath.dat','w');
    for ii = 1:N
        fprintf(fid, '%f\n', bath1D(ii));
    end
    fclose(fid);
end

%% Generate bathymetry edge
if iedge == 1
    edgex = bathEdge_iiP; edgey = bathEdge_jjP;
    dim = [size(edgex,1) size(edgex,2)];
    N = dim(1)*dim(2);
    % remove NaNs in edges
    for ii = 1:dim(1)-1
        for jj = 1:dim(2)-1
            if isnan(edgex(ii,jj))
                edgex(ii,jj) = max(data(ii,jj), data(ii+1,jj));
            end
            if isnan(edgey(ii,jj))
                edgey(ii,jj) = max(data(ii,jj), data(ii,jj+1));
            end
        end
    end
    for ii = 1:dim(1)
        if isnan(edgey(ii,end))
            edgey(ii,end) = data(ii,end);
        end
    end
    for ii = 1:dim(1)-1
        if isnan(edgex(ii,end))
            edgex(ii,end) = max(data(ii,end),data(ii+1,end));
        end
    end
    edgex(end,end) = data(end,end);
    for ii = 1:dim(2)
        if isnan(edgex(end,ii));
            edgex(end,ii) = data(end,ii);
        end
    end
    for ii = 1:dim(2)-1
        if isnan(edgey(end,ii));
            edgey(end,ii) = max(data(end,ii),data(end,ii+1));
        end
    end
    edgey(end,end) = data(end,end);
    aa = isnan(edgex);      bb = isnan(edgey);
    if sum(aa(:)) > 0 || sum(bb(:)) > 0
        error('NaNs in edges detected!');
    end
    edgex1D = reshape(edgex, N, 1);
    edgey1D = reshape(edgey, N, 1);
    % write x edge
    fid = fopen('edgex.dat','w');
    for ii = 1:N
        fprintf(fid, '%f\n', edgex1D(ii));
    end
    fclose(fid);
    % write y edge
    fid = fopen('edgey.dat','w');
    for ii = 1:N
        fprintf(fid, '%f\n', edgey1D(ii));
    end
    fclose(fid);
end

%% Generate tide file
if itide == 1
    load(fnameTide);
    N = size(tideout,1);
    fid = fopen('tide.dat','w');
    for ii = 1:N
        t = tNum + 86400*(tideout(ii,1) - tNumM);
        fprintf(fid, '%f %f\n', [round(t) tideout(ii,2)]);
    end
    fclose(fid);
end

%% Generate inflow file
if iinflow == 1
%     load('dailypump.mat');
    load(fnamePump);
    N = size(pumpout,1);
    fid = fopen('inflow.dat','w');
    for ii = 1:N
        t = tNum + 86400*(pumpout(ii,1) - tNumM);
        fprintf(fid, '%f %f\n', [round(t) pumpout(ii,2)]);
    end
    fclose(fid);
end

%% Generate wind file
if iwind == 1
    load('windspddir.mat');
    N = size(windoutV2,1);
    fid = fopen('windspd.dat','w');
    for ii = 1:N
        t = tNum + 86400*(windoutV2(ii,1) - tNumM);
        fprintf(fid, '%f %f\n', [round(t) windoutV2(ii,2)]);
    end
    fclose(fid);

    fid = fopen('winddir.dat','w');
    for ii = 1:N
        t = tNum + 86400*(windoutV2(ii,1) - tNumM);
        fprintf(fid, '%f %f\n', [round(t) windoutV2(ii,3)]);
    end
    fclose(fid);
end

%% Generate salinity file
if isalinity == 1
    load('salinity_BC.mat');
    N = size(salout,1);
    fid = fopen('salinityBC.dat','w');
    for ii = 1:N
        t = tNum + 86400*(salout(ii,1) - tNumM);
        fprintf(fid, '%f %f\n', [round(t) salout(ii,2)]);
    end
    fclose(fid);
end
%% Generate scalar ic
if iscalaric == 1
    load(fnameic);
%     scalar = Zhat;
    dim = [size(scalar,1) size(scalar,2)];
    N = dim(1)*dim(2);
    scalar1D = reshape(scalar, N, 1);
    fid = fopen('scalaric.dat','w');
    for ii = 1:N
        fprintf(fid, '%f\n', scalar1D(ii));
    end
    fclose(fid);
end

%% Generate subgrid data
if isubgrid == 1
    load(fnamesubA);
    load(fnamesubB);
    % get the number of elements to be written
    N = size(subB.bottom,1) * size(subB.bottom,2) * length(subA.surf);
    M = size(subB.bottom,1) * size(subB.bottom,2);
    subA = rmfield(subA, {'surf', 'dx', 'Dx'});
    if isfield(subA,{'block'})
        subA = rmfield(subA, {'block'});
    end
    fieldA = fieldnames(subA);
    fieldB = fieldnames(subB);
    % write to file
    mkdir('subdata');
    cd('subdata');
    for kk = 1:length(fieldB)
        field = fieldB{kk};
        bath1D = reshape(subB.(field), M, 1);
        fid = fopen([char(field),'.dat'], 'w');
        if strcmp(field, 'Exflag') || strcmp(field, 'Eyflag')
            for ii = 1:M
                fprintf(fid, '%d\n', bath1D(ii));
            end
        else
            for ii = 1:M
                fprintf(fid, '%f\n', bath1D(ii));
            end
        end
        fclose(fid);
    end
    for kk = 1:length(fieldA)
        field = fieldA{kk};
        area1D = reshape(subA.(field), N, 1);
        fid = fopen([char(field),'.dat'], 'w');
        for ii = 1:N
            fprintf(fid, '%f\n', area1D(ii));
        end
        fclose(fid);
    end
    cd ..;
end






