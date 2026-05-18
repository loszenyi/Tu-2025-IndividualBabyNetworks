clear
master_dir = '/data/wheelock/data1/people/Cindy/rfMRI_VAE/rsfMRI-VAE'
addpath(master_dir)
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup'));
rmpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup/tmp');
rmpath('/data/wheelock/data1/people/Cindy/PublicRepo/BrBx-HSB_infomap_cleanup/legacy');
addpath(genpath('/data/wheelock/data1/people/Cindy/PublicRepo/NLA_toolbox_070319/NLA visualizationfunctions'));
addpath([master_dir,'/subjlists']);
addpath([master_dir,'/CIFTI_read_save'])

load('/data/wheelock/data1/people/Cindy/PublicRepo/GeneralizabilityAreaParcellation/NetworkCommunityAssignments/IM_Tu326_12Networks.mat','IM');
IM_Tu12 = IM;
load('/data/wheelock/data1/people/Cindy/PublicRepo/GeneralizabilityAreaParcellation/NetworkCommunityAssignments/IM_Tu326_19Networks.mat','IM');
IM_Tu19 = IM;
load('IM_Gordon_13nets_333Parcels_Renamed2024.mat','IM');
IM_Gordon = IM;
load('/data/wheelock/data1/parcellations/IM/IM_MyersLabonte_23Networks.mat','IM');
IM_MyersLabonte = IM; IM_MyersLabonte.Nets{end} = 'Unspecified';
clear IM
% Load parcels
Parcels_Tu_326 = smartload('Parcels_Tu_326.mat');
Parcels_Gordon = smartload('Parcels_Gordon.mat');
Parcels_MyersLabonte283 = smartload('Parcels_Myers_Labonte_202310.mat');

load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;clear MNIl MNIr

parcelname = 'Tu_326'
load('eLABE_Y2_Y3_template_dim326_centroids.mat');

IM = IM_Tu19;
Parcels = Parcels_Tu_326;

clear C_bestk assn_bestk
C_bestks = who('C_bestk*')
bestks = reshape(cellfun(@str2double,regexp(C_bestks,'\d+','match')),1,[]);
Coptimized = single(Coptimized);
centroids_init_corr  = single(centroids_init_corr);
for bestk = bestks
    eval(['C_bestk',num2str(bestk),'= single(C_bestk',num2str(bestk),');'])
end
%% Assign each vertex
namestr = 'BCP_Jan2023_QCpass_asleep_atleast8min_UNC_UMN_20240124'%'eLABE_Y2_N113_atleast600frames'%'BCP_Jan2023_QCpass_awake_atleast8min_UNC_UMN_20240124'% 'eLABE_Y2_N113_atleast600frames'%'BCP_Jan2023_QCpass_asleep_atleast8min_UNC_UMN_20240124'
result_file_path = ['./template_matching_results/',namestr,'_wholesession.mat'];
T = readtable([master_dir,'/subject_tables/',namestr,'_vars.csv']);
YearGroup = T.ses_id{1}
to_overwrite = 0;

subs =  importdata([namestr,'.txt']);
if contains(namestr,'BCP')
    cohortfile = ['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/cohortfiles/cohortfiles_',namestr,'_2.55sigma.txt'];
    tmasklist=['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/tmasklist/tmasklist_',namestr,'_2.55sigma_QC_and_FDpt2_removeoutlierwholebrain_outliercalculatedonlowFDframes.txt'];
else
    cohortfile = ['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/cohortfiles/cohortfiles_',namestr,'.txt'];
    tmasklist=['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/tmasklist/tmasklist_',namestr,'.txt'];
end
% Read in subject names, functional volume locations, and surface directory
[subjects, cifti_files,~,~] = textread(cohortfile,'%s%s%s%s');
% Read in tmasks
[tmasksubjects, tmaskfiles]=textread(tmasklist,'%s%s');
assert(isequal(tmasksubjects,subjects),'tmasklist subjects do not match cohortfile subjects')
Nsubs = length(subs);
if exist(result_file_path,'file') 
    disp('Existing results loaded');
    load(result_file_path); % so we only write the new data
    if (to_overwrite==1)
        warning('Are you sure you want to overwrite the data?')
    end
else    
end
if ~exist(result_file_path,'file') || (to_overwrite ==1)
    tic
    for i = 1:Nsubs 
        fprintf('loading %i out of %i sessions...\n',i,Nsubs);
        subjectname = subs{i};
        tmask = importdata(tmaskfiles{i});
        
        if contains(T.dataset{1},'eLABE')
            [~,dtseriesname] = fileparts(cifti_files{i});
            dtseriescifti = ft_read_cifti_mod(['/data/wheelock/data1/datasets/eLABE/dtseries/imputedbyneighbors/',YearGroup,'/',dtseriesname,'.nii']);
        else
            dtseriescifti = ft_read_cifti_mod(cifti_files{i});
        end
        brainstructure = dtseriescifti.brainstructure(dtseriescifti.brainstructure>0);
        dtseries = single(dtseriescifti.data);
        clear dtseriescifti
        dtseries = dtseries(brainstructure<3,tmask==1);
        
        if contains(T.dataset{1},'BCP')
            ptseries_filename = dir(['/data/wheelock/data1/datasets/BCP/January2023/ptseries/',parcelname,'/',subjectname,'*.ptseries.nii']);
        elseif contains(T.dataset{1},'eLABE')
            ptseries_filename = dir(['/data/wheelock/data1/datasets/eLABE/ptseries/',YearGroup,'/',parcelname,'/',subjectname,'*.ptseries.nii']);
        end
        ptseries = ft_read_cifti_mod(fullfile(ptseries_filename.folder,ptseries_filename.name));
        ptseries = single(ptseries.data(:,tmask==1));
        
        data =  corr(dtseries',ptseries');
        
        for bestk = bestks
            eval(['C_bestk = C_bestk',num2str(bestk),';'])
            [D_top2,assn_top2] = pdist2(C_bestk,data,'correlation','Smallest',2);
            eval(['assn_all.bestk',num2str(bestk),'(:,i,:) = [assn_top2]'';']);
            eval(['D_all.bestk',num2str(bestk),'(:,i,:) = [D_top2]'';']);
        end

        [D_top2,assn_top2] = pdist2(centroids_init_corr,data,'correlation','Smallest',2); % Find the nearest centroid
        assn_all.init_corr(:,i,:) =  [assn_top2]';
        D_all.init_corr(:,i,:) = [D_top2]';
        toc
    end
    save(result_file_path,'assn_all*','D_all*');
end
return

