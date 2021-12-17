function step5_interpolateRegions(hemiMirror,parcellation_name)


%hemiMirror=
%0-Do not mirror samples between hemisphere
%1-Mirror samples and consider only left hemisphere
%2-Mirror samples between hemispheres (left samples are fliped and
%considered in the right hemisphere and viceversa)


        
path_Allen='AIBS_map/';

%Select mirroring method
if hemiMirror==0
    path_genes_general=[path_Allen 'gene_matrices/'];
elseif hemiMirror==1
    path_genes_general=[path_Allen 'gene_matrices_hemi_mirror/'];
elseif hemiMirror==2
    path_genes_general=[path_Allen 'gene_matrices_mirror/'];
end

[s1 s2 s3]=mkdir(path_genes_general);

interpMethod=2;
%Interpolation method, in case there are regions without
%samples 有些区域没有样本的情况下，采用插值法
if interpMethod==1
    interpMethodName='nearest';
elseif interpMethod==2
    interpMethodName='linear';
elseif interpMethod==3
    interpMethodName='extreme';
end


dir_img=dir([path_genes_general 'raw/']);
parc_cortex=[];gene_samples=[];samples_coor_mni=[];samples_coor_vox=[];

%Estimate which samples lies within each region for each donor
%每个donor在每个区域内有哪些样本
folder_name_list={'normalized_microarray_donor9861','normalized_microarray_donor10021','normalized_microarray_donor12876','normalized_microarray_donor14380','normalized_microarray_donor15496','normalized_microarray_donor15697'};
% folder_name_list={'normalized_microarray_donor9861','normalized_microarray_donor10021','normalized_microarray_donor12876','normalized_microarray_donor14380','normalized_microarray_donor15697'};
for img=1:numel(folder_name_list)
    donor_name=folder_name_list{img};
    %Voxels that include a sample
    path_subject=[path_genes_general 'raw/' donor_name '/'];
    samples_coor_vox_single=load([path_subject 'samples_coor_vox.mat']);
    samples_coor_vox_single=samples_coor_vox_single.samples_coor_vox;
    %Add one to each index because matlab is 1-based
    samples_coor_vox_single=samples_coor_vox_single+ones(size(samples_coor_vox_single));
    samples_coor_vox=[samples_coor_vox;samples_coor_vox_single];
    
    gene_samples_single=load([path_subject 'gene_samples.mat']);%gene_samples_single是结构体，里面存放gene_samples，gene_samples是1096*20647（对于9861被试）矩阵
    
    %mni coordinates of the samples, just in case need to be interpolated
    samples_coor_mni_single=load([path_subject 'samples_coor_mni.mat']);
    samples_coor_mni_single=samples_coor_mni_single.samples_coor_mni;
    samples_coor_mni=[samples_coor_mni;samples_coor_mni_single];
    
    numGenes=size(gene_samples_single.gene_samples,2);%返回矩阵列数，这里是返回20647
    gene_samples=[gene_samples;gene_samples_single.gene_samples];%gene_samples为1096*20647（对于被试9861）
    donor_name=donor_name(23:end);
    path_parcellation=[path_Allen 'Allen_FS/' donor_name '/parcellation/'];
    vol_parc=load_nifti([path_parcellation parcellation_name]);%结构体
    vol_samples=vol_parc;
    vol_parc=vol_parc.vol;%256*2556*256矩阵
    vol_samples.vol=zeros(size(vol_samples.vol));%生成全0矩阵
    vol_cortical=vol_parc;
    
    centroids_parcels=centroidFromParcellationFun([path_Allen 'Allen_FS/' donor_name '/parcellation/'], parcellation_name);%308*3矩阵，坐标

    %Total parcel depend whether we choose one or both hemisphere
    if hemiMirror==1
        npar=sum(centroids_parcels(:,1)<0); 
    else
        npar=size(centroids_parcels,1);%返回矩阵行数，308
    end
    npar_vol=numel(unique(vol_parc))-1;%unique(vol_parc)得到的是309*1的矩阵，内容从0-309
    vol_cortical(find(vol_cortical>npar))=0;
    
    %Create a nifty file showing the location of each sample in each donor
    %brain
    for ic=1:size(samples_coor_vox_single,1)
        parc_cortex(end+1)=vol_parc(samples_coor_vox_single(ic,1),samples_coor_vox_single(ic,2),samples_coor_vox_single(ic,3));
        value_par=vol_parc(samples_coor_vox_single(ic,1),samples_coor_vox_single(ic,2),samples_coor_vox_single(ic,3));
        for ind1=-1:1:1
            for ind2=-1:1:1
                for ind3=-1:1:1
                    vol_samples.vol(samples_coor_vox_single(ic,1)+ind1,samples_coor_vox_single(ic,2)+ind2,samples_coor_vox_single(ic,3)+ind3)=value_par;
                end
            end
        end
    end
    parc_cortex(find(parc_cortex>npar))=0;
    %Location of the samples for Quality control. It should overlap with
    %the file donors_name/mri/T1.nii
    save_nifti(vol_samples,[path_Allen 'samples_location_' donor_name '.nii.gz']);
    
end
%Calculate proportion of the samples used
if hemiMirror==2
    numSamples=numel(parc_cortex)/2;
    ratio1=parc_cortex(1:numSamples)>0;
    ratio2=parc_cortex(1+numSamples:end)>0;
    parc_cortex_ratio=sum((ratio1+ratio2)>0)/numSamples;
else
    parc_cortex_ratio=sum(parc_cortex>0)/numel(parc_cortex);
    numSamples=numel(parc_cortex);
end
parc_cortex_samples=zeros(npar,1);
gene_regional_expression=nan(npar,numGenes);
%Extract all the samples associated to each parcel
for ip=1:npar
    samplesRegion=find(parc_cortex==ip);%ip与parc_cortex值是否相等，相等返回1，不等返回0；find函数是返回矩阵或者向量中不是0的元素的位置索引。
    samplesRegion(find(samplesRegion>numSamples))=samplesRegion(find(samplesRegion>numSamples))-numSamples;
    if isempty(samplesRegion)
        if interpMethod==1
            gene_regional_expression(ip,:)=allen_interp_nearest(gene_samples,samples_coor_mni,centroids_parcels(ip,:));
        elseif interpMethod==2
            gene_regional_expression(ip,:)=allen_interp_linear_interp(gene_samples,samples_coor_mni,centroids_parcels(ip,:));
        elseif interpMethod==3
            gene_regional_expression(ip,:)=9999*ones(1,size(gene_regional_expression,2));
        else
            error
        end
    else
        %Calculate the gene expression of a region as the median gene
        %expression values
        gene_regional_expression(ip,:)=median(gene_samples(samplesRegion,:),1);
        parc_cortex_samples(ip)=numel(samplesRegion);
    end
end

%Save files
path_output=path_genes_general;


[s1 s2 s3]=mkdir(path_output);
save([path_output 'gene_regional_expression.mat'],'gene_regional_expression');
gene_regional_expression=zscore(gene_regional_expression);
save([path_output 'gene_regional_expression_zscored.mat'],'gene_regional_expression');
gene_regional_correlations=corrcoef(gene_regional_expression');
save([path_output 'gene_regional_correlations.mat'],'gene_regional_correlations');
save([path_output 'parc_cortex_samples.mat'],'parc_cortex_samples');
save([path_output 'parc_cortex_ratio.mat'],'parc_cortex_ratio');
display(['Expression matrices created:' path_output 'gene_regional_expression.mat']);
