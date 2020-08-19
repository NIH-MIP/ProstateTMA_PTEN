function[] = spatial_process_training(xmlLoc,saveFolder)

        %csv with list of tmas by id with xml file name (if exists)
        all_info = readtable('/path/to/xml_list.txt');
        
        tma_info = readtable('/path/to/allTMAs_info.txt');
        tma_info = table2cell(tma_info);

        for i = 1:size(tma_info,1)
            tma_info{i,9} = ['TMA' int2str(tma_info{i,1}) '_' int2str(tma_info{i,2}) '_' int2str(tma_info{i,3}) '_' tma_info{i,5}];
        end

        for i = 1:size(all_info,1)
            tma_id = all_info{i};
            if(~exist([saveFolder filesep tma_id]))
                mkdir([saveFolder filesep tma_id])
            end

            if(~exist([saveFolder filesep tma_id filesep 'spatial']))
                mkdir([saveFolder filesep tma_id filesep 'spatial'])
                mkdir([saveFolder filesep tma_id filesep 'spatial' filesep 'loss'])
                mkdir([saveFolder filesep tma_id filesep 'spatial' filesep 'intact'])
            end


            case_info = tma_info(strcmpi(tma_info(:,9),tma_id),:);
            caseData = bfGetReader([mainDir filesep case_info{1,6}]);
            caseMeta = caseData.getMetadataStore();
            caseData.setSeries(0);
            imgsizex = eval(caseMeta.getPixelsSizeX(0));
            imgsizey = eval(caseMeta.getPixelsSizeY(0));
            if((case_info{1,7}+tmaSz)>imgsizex)
                imgsizex = imgsizex-case_info{1,7};
            else
                imgsizex=tmaSz;
            end
            if((case_info{1,8}+tmaSz)>imgsizey)
                imgsizey = imgsizey-case_info{1,8};
            else
                imgsizey = tmaSz;
            end
            I1 = bfGetPlane(caseData,1,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
            I2 = bfGetPlane(caseData,2,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
            I3 = bfGetPlane(caseData,3,case_info{1,7}+1,case_info{1,8}+1,imgsizex,imgsizey);
            tma_img(:,:,1) = I1;
            tma_img(:,:,2) = I2;
            tma_img(:,:,3) = I3;
            
            if(strcmpi(all_info{i,2},'None'))
                spatial_noxml(tma_img, tma_id, saveFolder)
            else
                spatial_hasxml(tma_img,tma_id,all_info{i,2},xmlLoc,saveFolder)
            end
            clear I1 I2 I3 tma_img
        end
end




function[] = spatial_hasxml(tma_img,tma_id,tma_xml_file,xmlLoc,saveFolder)

         tma_info = readtable('/path/to/tma_info.txt');
         tma_info = table2cell(tma_info);
         boxSz = 50;
%          saveFolder = '/path/to/save';
         for i = 1:size(tma_info,1)
             tma_info{i,9} = ['TMA' tma_info{i,1} '_' int2str(tma_info{i,2}) '_' int2str(tma_info{i,3}) ];
         end

        tma_outcome = tma_info{find(strcmpi(tma_info(:,6), strrep(tma_xml_file,'xml','svs'))),5};
%         tma_id = [strrep(tma_xml_file,'.xml','')];     
        xml_file = [xmlLoc filesep tma_xml_file];
        xy = xml_parse(xml_file);
        
        null_mask = zeros([size(tma_img,1) size(tma_img,2)]);
        for region = 1:size(xy,2)
            pts = xy{1,region};
            if(pts(1,1) == 0 && pts(1,2) == 0)
            elseif(pts(end,1) == 0 && pts(end,2) == 0)
            else
                rmask = poly2mask(pts(:,1),pts(:,2),size(tma_img,1),size(tma_img,2));
                rinds = find(rmask>0);
                null_mask(rinds) = 1;
            end    
        end
        

        bw_img = double(rgb2gray(tma_img))./255;
        rgb = imresize(tma_img,[200 200]);
        gray_image = rgb2gray(rgb);
        bw = imbinarize(gray_image);
        mask = zeros(size(gray_image));
        mask(find(bw==0))=1;
        bw = imfill(mask,'holes');
        se = strel('disk',10);
        bw = imclose(bw,se);
        se = strel('disk',5);
        bw = imdilate(bw,se);
        % imshow(bw)
        tissue = imresize(bw,size(bw_img),'nearest');


        inds_null = find(null_mask>0);
        inds_tissue = find(tissue>0);
        inds_pos = setdiff(inds_tissue,inds_null);
        pos_mask = zeros(size(bw_img));
        pos_mask(inds_pos) = 1;
        
        null_rois = bwlabel(null_mask);
        numrois = unique(null_rois(null_rois>0));
        
        %go through simple tissue mask first, only take 50x50 boxes that
        %are either 90% PTEN loss or 90% in PTEN normal tissue
        
        % start with tissue
        for j = 1:numel(numrois)
            jinds = find(null_rois==numrois(j));
            [jroi(:,1),jroi(:,2)]= ind2sub(size(bw_img),jinds);
            y_min = round(min(jroi(:,2))); y_max = round(max(jroi(:,2)));
            x_min = round(min(jroi(:,1))); x_max = round(max(jroi(:,1)));
            [xmin, xmax, xstride] = findminmax(x_min, x_max, size(bw_img,1), boxSz);
            [ymin, ymax, ystride] = findminmax(y_min, y_max, size(bw_img,2), boxSz);
            
            %jmask = poly2mask(jroi(:,2),jroi(:,1),size(bw_img,1),size(bw_img,2));
            jmask = zeros(size(null_rois));
            jmask(jinds) = 1;
            x_samp = ceil((xmax-xmin)/boxSz);
            y_samp = ceil((ymax-ymin)/boxSz);
            
            box_count = 1;
            for xi = 1:x_samp
                for yi = 1:y_samp
                    x_loc = xmin + (xi-1)*boxSz - (xi-1)*ceil(xstride);
                    y_loc = ymin + (yi-1)*boxSz - (yi-1)*ceil(ystride);
                    sm_mask = jmask(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1);
                    box_mask = zeros(size(jmask));
                    box_mask(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1) = 1;
                    perc_in = length(find(sm_mask>0))/(boxSz*boxSz);
                    if(perc_in >= 0.9)
                        bw_sm  = bw_img(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                        find_white = find(bw_sm > 0.8);
                        prop_white = length(find_white)/(size(bw_sm,1)*size(bw_sm,2));
                        if(prop_white < 0.95)
                            %return
                            tma_crop = tma_img(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                            imwrite(tma_crop,[saveFolder filesep tma_id filesep 'spatial' filesep 'loss' filesep strrep(all_tmas{i},'.xml','') '_roi' int2str(j) '_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' tma_outcome '.png']);
                            disp([tma_id ' box ' int2str(box_count) ' ' int2str(x_loc) ' ' int2str(y_loc)])
                            box_count = box_count + 1;
                        end
                    end
                end
            end
            clear jroi
        end
        
        % now for normal
        [subs_pos(:,1),subs_pos(:,2)]  = ind2sub(size(bw_img),inds_pos);
        y_min = round(min(subs_pos(:,2))); y_max = round(max(subs_pos(:,2)));
        x_min = round(min(subs_pos(:,1))); x_max = round(max(subs_pos(:,1)));
        [xmin, xmax, xstride] = findminmax(x_min, x_max, size(bw_img,1), boxSz);
        [ymin, ymax, ystride] = findminmax(y_min, y_max, size(bw_img,2), boxSz);
        x_samp = ceil((xmax-xmin)/boxSz);
        y_samp = ceil((ymax-ymin)/boxSz);
        box_count = 1;
        for xi = 1:x_samp
            for yi = 1:y_samp
                x_loc = xmin + (xi-1)*boxSz - (xi-1)*ceil(xstride);
                y_loc = ymin + (yi-1)*boxSz - (yi-1)*ceil(ystride);
                sm_mask = pos_mask(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1);
                perc_in = length(find(sm_mask>0))/(boxSz*boxSz);
                if(perc_in >= 0.9)
                    bw_sm  = bw_img(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                    find_white = find(bw_sm > 0.8);
                    prop_white = length(find_white)/(size(bw_sm,1)*size(bw_sm,2));
                    if(prop_white < 0.9)
                        %return
                        tma_crop = tma_img(x_loc:x_loc+boxSz-1,y_loc:y_loc+boxSz-1,:);
                        imwrite(tma_crop,[saveFolder filesep tma_id filesep 'pos' filesep strrep(all_tmas{i},'.xml','') '_tissue_box' int2str(box_count) '-' int2str(x_loc) '-' int2str(y_loc) '_' tma_outcome '.png']);
                        disp([tma_id ' box ' int2str(box_count) ' ' int2str(x_loc) ' ' int2str(y_loc)])
                        box_count = box_count + 1;
                    end
                end
            end
        end
        clear subs_pos tma_img
%     end
end

function [] = spatial_noxml(tma_img, tma_id, saveFolder)

    saveDir = [saveFolder filesep tma_id filesep 'spatial' filesep 'intact'];

    imgSize = 50;
    tmaSize = [size(tma_img,1) size(tma_img,2)];
    stride = 50;
    
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
            if(wscount<0.9)
                %disp([tma_id '   sub: ' int2str(sub_count)])
                imwrite(sub,[saveDir filesep tma_id '_sub' int2str(sub_count) '.png']);
                sub_inds = cat(1, sub_inds, {sub_count} ,{new_sub});
                sub_count = sub_count+1;
            end
            clear sub bwsub wscount
        end
    end
     
    %save([saveSub filesep tma_id '_linkers.mat'],'sub_inds','-v7.3');  
    
end

function[xy] = xml_parse(xml_file)
   
    xDoc = xmlread(xml_file); %get xml doc
    Regions=xDoc.getElementsByTagName('Region'); % get a list of all the region tags

    %get region labels
    for regioni = 0:Regions.getLength-1 %region tags start at 0
        Region=Regions.item(regioni);  %for each region tag
        label1 = Region.getAttribute('NegativeROA'); %ground truth label;
        label2 = Region.getAttribute('Text');
        verticies=Region.getElementsByTagName('Vertex'); %get a list of all the vertexes (which are in order)
        xy{1,regioni+1}=zeros(verticies.getLength-1,2); %allocate space for them
        xy{2,regioni+1}=char(label1); %grab label for verticies 
        xy{3,regioni+1}=char(label2);
        for vertexi = 0:verticies.getLength-1 %iterate through all verticies
            x=str2double(verticies.item(vertexi).getAttribute('X')); %get the x value of that vertex
            y=str2double(verticies.item(vertexi).getAttribute('Y')); %get the y value of that vertex
            xy{1,regioni+1}(vertexi+1,:)=[x,y]; % finally save them into the array
        end   
    end

    
end

function [minout, maxout, stride] = findminmax(minval, maxval, sizeval, boxSz)

    leftover = round((ceil((maxval-minval)/boxSz)-((maxval-minval)/boxSz))*boxSz); 
    stride = 0;

    if(mod(leftover,2) == 1) %odd
        % if there is space to the left for leftover/2
        if((minval - ceil(leftover/2)) > 0)
            % if there is space to the right for leftover/2
            if((maxval + ceil(leftover/2)) < sizeval)
                minval = minval - ceil(leftover/2);
                maxval = maxval + ceil(leftover/2)-1;
            else 
                %otherwise go
                minval = minval - (leftover - (sizeval - maxval));
                maxval = sizeval;
                if(minval < 1)
                    minval = 1;
                    leftover = round((ceil((maxval-minval)/boxSz)-((maxval-minval)/boxSz))*boxSz)+1;
                    stride = leftover/(ceil((maxval-minval)/boxSz)-1);
                    %stride = mod((maxval-minval +1),boxSz)/(floor((maxval-minval+1)/boxSz)-1);
                end
            end
        else
            maxval = maxval + (leftover - minval + 1);
            minval = 1;
            %if the other side is also out of space
            if(maxval > sizeval)
                maxval = sizeval;
                leftover = round((ceil((maxval-minval)/boxSz)-((maxval-minval)/boxSz))*boxSz)+1;
                stride = leftover/(ceil((maxval-minval)/boxSz)-1);
                %stride = mod((maxval-minval+1),boxSz)/(floor((maxval-minval+1)/boxSz)-1);
            end 
        end
    else %even 
        %if min val has enough space to move leftover/2
        if((minval - leftover/2) > 0)
            %if maxval has enough space to move leftover/2
            if((maxval + leftover/2) < sizeval)
                minval = minval - leftover/2;
                maxval = maxval + leftover/2;
            %if not, set the difference to the other side
            else
                minval = minval - (leftover - (sizeval - maxval));
                maxval = sizeval;
                %if the other side is also out of space
                if(minval < 1)
                    minval = 1;
                    leftover = round((ceil((maxval-minval)/boxSz)-((maxval-minval)/boxSz))*boxSz)+1;
                    stride = leftover/(ceil((maxval-minval)/boxSz)-1);
                    %stride = mod((maxval-minval+1),boxSz)/(floor((maxval-minval+1)/boxSz)-1);
                end
            end
        %if the min side doesnt have enough, set the difference to the other side    
        else
            maxval = maxval + maxval + (leftover - minval + 1);
            minval = 1;
            %if the other side is also out of space
            if(maxval > sizeval)
                maxval = sizeval;
                leftover = round((ceil((maxval-minval)/boxSz)-((maxval-minval)/boxSz))*boxSz)+1;
                stride = leftover/(ceil((maxval-minval)/boxSz)-1);
                %stride = mod((maxval-minval+1),boxSz)/(floor((maxval-minval+1)/boxSz)-1);
            end            
        end
    end

    minout = minval;
    maxout = maxval;
    %disp(num2str(mod((maxout-minout),boxSz)))
end

