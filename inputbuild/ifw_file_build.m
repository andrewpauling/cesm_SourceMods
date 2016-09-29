%ifw_file_build

close all
clear all

%Specify the type of run. Can be:
    % Summer 
    % Winter
    % mthavg
% where mthavg is for constant or ramped freshwater forcing

runtype = 'summer';

% Set output file name
outname = ['ifw2000Gt_',runtype,'_FRESHWATER_gx1v6.nc'];

snc_setup; 
disp('--------------------------------------------')
disp('--------BEGINNING INPUT FILE BUILD----------')
disp('--------------------------------------------')
disp('                                            ')
disp('Reading CESM variables...')

load cesm_vars.mat

disp('Done reading')
disp(' ')
disp('Loading mask data...')

load('re_max_level_gx1v6.mat');
load('refixed_layers_gx1v6.mat');
load('refixed_coast_mask_gx1v6.mat');
zero_ind = find(isnan(final_max)==1);
other_ind = find(final_max == -10);
final_max(other_ind) = 0;
final_max(zero_ind) = 0;
disp('Done loading')
disp(' ')

cell_area3d = zeros(1,384,320);
cell_area3d(1,:,:) = cell_area;
sec_yr = 365*24*3600;
Gt_to_kg = 1e-12;

%Create an empty netCDF file into which the data will be put.
outfile = outname;

disp(['Creating empty output file: ',outname,' ...'])
mode = bitor(nc_clobber_mode,nc_64bit_offset_mode);
nc_create_empty(outfile,mode);
nc_adddim(outfile,'Z',60);
nc_adddim(outfile,'Y',384);
nc_adddim(outfile,'X',320);
disp('Done')
disp(' ')
disp('Beginning hosing variable build...')

if strcmp(runtype,'mthavg')==1;
    
    melt_total = 2000;
    
    for nn = 1:12;
        FRESHWATER = ifw_fill(melt_total,nn,new_coast,my_layers);
        %Use function ifw_fill to create 60x384x320 variable containing
        %cells into which freshwater will be put and the amount for the
        %given month in Gt/month.
        for kk = 1:60;
            FRESHWATER(kk,:,:) = FRESHWATER(kk,:,:)./cell_area3d./Gt_to_kg;
        end
        %Divide by cell area and convert from Gt to kg.
        
        FRESHWATER = FRESHWATER.*conversion;
        FRESHWATER = FRESHWATER./sec_yr;
        %convert to g/g/yr from kg/m^2/yr then to g/g/s.
        
        varstruct.Name = ['FRESHWATER_',num2str(nn)];
        varstruct.Datatype = 'float';
        varstruct.Attribute(1).Name = 'long_name';
        varstruct.Attribute(2).Name = 'units';
        varstruct.Attribute(3).Name = 'time_op';
        varstruct.Attribute(4).Name = 'coordinates';
        varstruct.Attribute(5).Name = 'missing_value';
        varstruct.Attribute(1).Value = ['Hosing Freshwater Flux in model'...
            ,' salt units per second'];
        varstruct.Attribute(2).Value = '(g/g)*(cm/s)';
        varstruct.Attribute(3).Value = 'average';
        varstruct.Attribute(4).Value = 'TLONG TLAT z_t';
        varstruct.Attribute(5).Value = -1e34;
        
        varstruct.Dimension = {'Z','Y','X'};
        nc_addvar(outfile,varstruct);
        start = [0 0 0]; count = [60,384 320];
        nc_varput(outfile,['FRESHWATER_',num2str(nn)],FRESHWATER,start,...
            count);
        
        disp(['Finished month',num2str(nn)]);
    end

elseif strcmp(runtype,'summer')==1;
    
    melt_total = 4000;
    
    for nn = [1:5,12];
        FRESHWATER = ifw_fill(melt_total,nn,new_coast,my_layers);
        %Use function ifw_fill to create 60x384x320 variable containing
        %cells into which freshwater will be put and the amount for the
        %given month in Gt/month.
        for kk = 1:60;
            FRESHWATER(kk,:,:) = FRESHWATER(kk,:,:)./cell_area3d./Gt_to_kg;
        end
        %Divide by cell area and convert from Gt to kg.
        
        FRESHWATER = FRESHWATER.*conversion;
        FRESHWATER = FRESHWATER./sec_yr;
        %convert to g/g/yr from kg/m^2/yr then to g/g/s.
        
        varstruct.Name = ['FRESHWATER_',num2str(nn)];
        varstruct.Datatype = 'float';
        varstruct.Attribute(1).Name = 'long_name';
        varstruct.Attribute(2).Name = 'units';
        varstruct.Attribute(3).Name = 'time_op';
        varstruct.Attribute(4).Name = 'coordinates';
        varstruct.Attribute(5).Name = 'missing_value';
        varstruct.Attribute(1).Value = ['Hosing Freshwater Flux in model'...
            ,' salt units per second'];
        varstruct.Attribute(2).Value = '(g/g)*(cm/s)';
        varstruct.Attribute(3).Value = 'average';
        varstruct.Attribute(4).Value = 'TLONG TLAT z_t';
        varstruct.Attribute(5).Value = -1e34;
        
        varstruct.Dimension = {'Z','Y','X'};
        nc_addvar(outfile,varstruct);
        start = [0 0 0]; count = [60,384 320];
        nc_varput(outfile,['FRESHWATER_',num2str(nn)],FRESHWATER,start,...
            count);
        
        disp(['Finished month ',num2str(nn)])
    end

    for nn = 6:11;
        FRESHWATER = zeros(60,384,320);
        
        varstruct.Name = ['FRESHWATER_',num2str(nn)];
        varstruct.Datatype = 'float';
        varstruct.Attribute(1).Name = 'long_name';
        varstruct.Attribute(2).Name = 'units';
        varstruct.Attribute(3).Name = 'time_op';
        varstruct.Attribute(4).Name = 'coordinates';
        varstruct.Attribute(5).Name = 'missing_value';
        varstruct.Attribute(1).Value = ['Hosing Freshwater Flux in model'...
            ,' salt units per second'];
        varstruct.Attribute(2).Value = '(g/g)*(cm/s)';
        varstruct.Attribute(3).Value = 'average';
        varstruct.Attribute(4).Value = 'TLONG TLAT z_t';
        varstruct.Attribute(5).Value = -1e34;
        
        varstruct.Dimension = {'Z','Y','X'};
        nc_addvar(outfile,varstruct);
        start = [0 0 0]; count = [60,384 320];
        nc_varput(outfile,['FRESHWATER_',num2str(nn)],FRESHWATER,start,...
            count);
        
        disp(['Finished month ',num2str(nn)])
    end

else
    
    melt_total = 4000;
    
    for nn = 6:11;
        FRESHWATER = ifw_fill(melt_total,nn,new_coast,my_layers);
        %Use function ifw_fill to create 60x384x320 variable containing
        %cells into which freshwater will be put and the amount for the
        %given month in Gt/month.
        for kk = 1:60;
            FRESHWATER(kk,:,:) = FRESHWATER(kk,:,:)./cell_area3d./Gt_to_kg;
        end
        %Divide by cell area and convert from Gt to kg.
        
        FRESHWATER = FRESHWATER.*conversion;
        FRESHWATER = FRESHWATER./sec_yr;
        %convert to g/g/yr from kg/m^2/yr then to g/g/s.
        
        varstruct.Name = ['FRESHWATER_',num2str(nn)];
        varstruct.Datatype = 'float';
        varstruct.Attribute(1).Name = 'long_name';
        varstruct.Attribute(2).Name = 'units';
        varstruct.Attribute(3).Name = 'time_op';
        varstruct.Attribute(4).Name = 'coordinates';
        varstruct.Attribute(5).Name = 'missing_value';
        varstruct.Attribute(1).Value = ['Hosing Freshwater Flux in model'...
            ,' salt units per second'];
        varstruct.Attribute(2).Value = '(g/g)*(cm/s)';
        varstruct.Attribute(3).Value = 'average';
        varstruct.Attribute(4).Value = 'TLONG TLAT z_t';
        varstruct.Attribute(5).Value = -1e34;
        
        varstruct.Dimension = {'Z','Y','X'};
        nc_addvar(outfile,varstruct);
        start = [0 0 0]; count = [60,384 320];
        nc_varput(outfile,['FRESHWATER_',num2str(nn)],FRESHWATER,start,...
            count);
        
        disp(['Finished month ',num2str(nn)])
    end

    for nn = [1:5,12];
        FRESHWATER = zeros(60,384,320);
        
        varstruct.Name = ['FRESHWATER_',num2str(nn)];
        varstruct.Datatype = 'float';
        varstruct.Attribute(1).Name = 'long_name';
        varstruct.Attribute(2).Name = 'units';
        varstruct.Attribute(3).Name = 'time_op';
        varstruct.Attribute(4).Name = 'coordinates';
        varstruct.Attribute(5).Name = 'missing_value';
        varstruct.Attribute(1).Value = ['Hosing Freshwater Flux in model'...
            ,' salt units per second'];
        varstruct.Attribute(2).Value = '(g/g)*(cm/s)';
        varstruct.Attribute(3).Value = 'average';
        varstruct.Attribute(4).Value = 'TLONG TLAT z_t';
        varstruct.Attribute(5).Value = -1e34;
        
        varstruct.Dimension = {'Z','Y','X'};
        nc_addvar(outfile,varstruct);
        start = [0 0 0]; count = [60,384 320];
        nc_varput(outfile,['FRESHWATER_',num2str(nn)],FRESHWATER,start,...
            count);
        
        disp(['Finished month ',num2str(nn)])
    end
    
end
    
disp(' ')
disp('Adding CESM variables...')
%Add T-grid areas
varstruct.Name = 'TAREA';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(3).Name = 'coordinates';
varstruct.Attribute(4).Name = '_FillValue';
varstruct.Attribute(5).Name = 'missing_value';
varstruct.Attribute(1).Value = 'area of T cells';
varstruct.Attribute(2).Value = 'centimeter^2';
varstruct.Attribute(3).Value = 'TLONG TLAT';
varstruct.Attribute(4).Value = 9.96920996838687e+36;
varstruct.Attribute(5).Value = 9.96920996838687e+36;

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'TAREA',cell_area,start,count);
disp('Added TAREA')

%Add U-grid areas
varstruct.Name = 'UAREA';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(3).Name = 'coordinates';
varstruct.Attribute(4).Name = '_FillValue';
varstruct.Attribute(5).Name = 'missing_value';
varstruct.Attribute(1).Value = 'area of U cells';
varstruct.Attribute(2).Value = 'centimeter^2';
varstruct.Attribute(3).Value = 'TLONG TLAT';
varstruct.Attribute(4).Value = 9.96920996838687e+36;
varstruct.Attribute(5).Value = 9.96920996838687e+36;

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'UAREA',uarea,start,count);
disp('Added UAREA')
clear varstruct

%Add T-grid latitudes
varstruct.Name = 'TLAT';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(1).Value = 'Latitude (T grid)';
varstruct.Attribute(2).Value = 'degrees north';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'TLAT',tlats,start,count);
disp('Added TLAT')

%Add T-grid longitudes
varstruct.Name = 'TLONG';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(1).Value = 'Longitude (T grid)';
varstruct.Attribute(2).Value = 'degrees east';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'TLONG',tlons,start,count);
disp('Added TLONG')

%Add U-grid latitudes
varstruct.Name = 'ULAT';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(1).Value = 'Latitude (U grid)';
varstruct.Attribute(2).Value = 'degrees north';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'ULAT',ulats,start,count);
disp('Added ULAT')

%Add U-grid longitudes
varstruct.Name = 'ULONG';
varstruct.Datatype = 'double';
varstruct.Attribute(1).Name = 'long_name';
varstruct.Attribute(2).Name = 'units';
varstruct.Attribute(1).Value = 'Longitude (U grid)';
varstruct.Attribute(2).Value = 'degrees east';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'ULONG',ulons,start,count);
disp('Added ULONG')

%Add region mask

varstruct.Name = 'REGION_MASK';
varstruct.Datatype = 'int';
varstruct.Attribute(1).Name = 'units';
varstruct.Attribute(2).Name = 'long_name';
varstruct.Attribute(3).Name = 'coordinates';
varstruct.Attribute(1).Value = 'Basin Index';
varstruct.Attribute(2).Value = 'basin index number (signed integers)';
varstruct.Attribute(3).Value = 'TLONG TLAT';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'REGION_MASK',mask,start,count);
disp('Added REGION_MASK')

%Add FWF_MAX_LEVEL

varstruct.Name = 'FWF_MAX_LEVEL';
varstruct.Datatype = 'int';
varstruct.Attribute(1).Name = 'units';
varstruct.Attribute(2).Name = 'long_name';
varstruct.Attribute(3).Name = 'coordinates';
varstruct.Attribute(4).Name = 'time_op';
varstruct.Attribute(1).Value = 'unitless';
varstruct.Attribute(2).Value = 'Maximum depth level for freshwater flux';
varstruct.Attribute(3).Value = 'TLONG TLAT';
varstruct.Attribute(4).Value = 'average';

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'FWF_MAX_LEVEL',final_max,start,count);
disp('Added FWF_MAX_LEVEL')

%Add KMT

varstruct.Name = 'KMT';
varstruct.Datatype = 'int';
varstruct.Attribute(1).Name = 'units';
varstruct.Attribute(2).Name = 'long_name';
varstruct.Attribute(3).Name = 'coordinates';
varstruct.Attribute(4).Name = '_FillValue';
varstruct.Attribute(5).Name = 'missing_value';
varstruct.Attribute(1).Value = 'unitless';
varstruct.Attribute(2).Value = 'k Index of Deepest Grid Cell on T Grid';
varstruct.Attribute(3).Value = 'TLONG TLAT';
varstruct.Attribute(4).Value = -2147483647;
varstruct.Attribute(5).Value = -2147483647;

varstruct.Dimension = {'Y','X'};
nc_addvar(outfile,varstruct);
start = [0 0]; count = [384 320];
nc_varput(outfile,'KMT',kmt,start,count);
disp('Added KMT')
disp(' ')
disp(['File saved at: ',outfile])
disp(' ')
disp('--------------------------------------------')
disp('---------INPUT FILE BUILD COMPLETE----------')
disp('--------------------------------------------')

