mainDir = '/data/MIP/PTEN/raw/external/raw';
%find tma outcome/class based on 
tma_info = readtable('/path/to/tma_info.txt');
tma_info = table2cell(tma_info);
for i = 1:size(tma_info,1)
    tma_info{i,9} = ['TMA' int2str(tma_info{i,1}) '_' int2str(tma_info{i,2}) '_' int2str(tma_info{i,3}) '_' tma_info{i,5}];
end
all_cases = tma_info(:,9);
tmaSz = 2000;
for ij = 1:size(all_cases,1)
%         try
        case_info = tma_info(strcmpi(tma_info(:,9),all_cases{ij}),:);
        caseData = bfGetReader([mainDir filesep case_info{1,6}]);
        caseMeta = caseData.getMetadataStore();
        caseData.setSeries(0);
        imgsizex = eval(caseMeta.getPixelsSizeX(0));
        imgsizey = eval(caseMeta.getPixelsSizeY(0));
        if(imgsizex > tmaSz )
            imgsizex = tmaSz ;
        end
        if(imgsizey > tmaSz )
            imgsizey = tmaSz ;
        end
        I1 = bfGetPlane(caseData,1,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
        I2 = bfGetPlane(caseData,2,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
        I3 = bfGetPlane(caseData,3,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
        tma_img(:,:,1) = I1;
        tma_img(:,:,2) = I2;
        tma_img(:,:,3) = I3;
        tma_id = case_info{1,9};
        disp(tma_id)
        %save two images, one large tma at full resolution and all subsequent boxes from smallers
        pullPatches_tma(tma_img, tma_id) 
%         catch ME
%             disp(['failed: ' batch_id ' ' int2str(bigcirc(ij,1)) '-' int2str(bigcirc(ij,2))])
%         end
clear I1 I2 I3 tma_img
end


function [] = pullPatches_tma(tma_img, tma_id)
    saveDir = ['/path/to/patches/dir/' filesep tma_id];
    saveSub = '/path/to/patches/dir/linkers';
    mkdir(saveDir)
    
    imgSize = 50;
    tmaSize = [size(tma_img,1) size(tma_img,2)];
    stride = 15;
    
    %find number we can fill space
    num_subs_updown = floor((tmaSize(1)-imgSize)/stride)+1;
    num_subs_leftright = floor((tmaSize(2)-imgSize)/stride)+1;

    sub_inds = [];
    sub_count = 1;
    for ii=1:num_subs_updown
        for jj = 1:num_subs_leftright
            new_sub = [stride*(ii-1)+1, stride*(ii-1)+imgSize, stride*(jj-1)+1, stride*(jj-1)+imgSize,];
            sub = tma_img(new_sub(1):new_sub(2),new_sub(3):new_sub(4),:);
            bwsub = double(rgb2gray(sub))./255;
            wscount = length(find(bwsub>0.8))./(size(bwsub,1)*size(bwsub,2));
            if(wscount<0.95)
                disp([tma_id '   sub: ' int2str(sub_count)])
                imwrite(sub,[saveDir filesep tma_id '_sub' int2str(sub_count) '.png']);
                sub_inds = cat(1, sub_inds, {sub_count} ,{new_sub});
                sub_count = sub_count+1;
            end
            clear sub bwsub wscount
        end
    end
     
    save([saveSub filesep tma_id '_linkers.mat'],'sub_inds','-v7.3');  
    
end
