import numpy as np
import openslide
import sys
import os
from PIL import Image
from xml.dom import minidom
import pandas as pd
from skimage import draw
import numpy as np
import matplotlib.pyplot as plt

class TMA_processFile:

    def __init__(self):
        self.csvname = '/path/to/info/file/tma_external_info.txt' #this file contains tma location and outcome as csv 
        self.data_location = '/location/of/raw/data' 
        self.save_location = '/location/of/saving/data'
        self.tma_size = (2000,2000) #size of tma core (this is avg size)
        self.base_size = 100


    def organize_tmas(self):
        df_csv = pd.read_csv(self.csvname,sep='\t',header=0)
        #print(df_csv['source'].dtypes)
        for index, tma in df_csv.iterrows():
            tma_id = "TMA" + str(tma['TMA']) + '_' + str(tma['row']) + '_' + str(tma['col'])
            tma_outcome = tma['PTEN']
            slide_name = tma['source']
            print(slide_name)
            startx = tma['xcoord']
            starty = tma['ycoord']
            if "Pronto" in slide_name:
                doSpatial = 0
            # else:
                # doSpatial = 1

            if os.path.isfile(os.path.join(self.data_location,slide_name)):
                self.read_tmaimage(slide_name=slide_name,startx=startx,starty=starty,tma_id=tma_id,tma_outcome=tma_outcome)
            if (doSpatial == 1):
                self.tma_doSpatial(slide_name=slide_name,startx=startx,starty=starty,tma_id=tma_id,tma_outcome=tma_outcome)


    def read_tmaimage(self,slide_name,startx,starty,tma_id, tma_outcome):
        oslide = openslide.OpenSlide(os.path.join(self.data_location,slide_name))
        shape = self.tma_size #width = oslide.dimensions[0]; #height = oslide.dimensions[1];
        patch = oslide.read_region((startx, starty), 0, shape);
        dir_id = tma_id + "_" +tma_outcome
        fname = os.path.join(self.save_location,'TMA',dir_id);
        os.mkdir(fname)
        print(fname)
        #patch.save(fname);

        #send to patches at varying levels
        os.mkdir(os.path.join(fname, 'TMA'))
        os.mkdir(os.path.join(fname,'5x'))
        os.mkdir(os.path.join(fname, '10x'))
        os.mkdir(os.path.join(fname, '20x'))
        if (doSpatial == 1):
            os.mkdir(os.path.join(fname, 'spatial'))

        #NOTE TO SELF THIS IS STUPID AND YOU SHOULD MAKE THESE VARIABLES IN INIT/CALCULATED
        sz5x = 4*self.base_size
        sz10x = 2*self.base_size
        boxSz = self.base_size
        patchnum = 1

        for x in range(1,int(img_size[0]/sz5x + 1)):
            for y in range(1,int(img_size[1]/sz5x + 1)):
                subpatch = patch.crop(box=(sz5x*(x-1),sz5x*(y-1),sz5x*x,sz5x*y))
                subpatch_name = tma_id + "_sub" + str(patchnum) + "_" + tma_outcome + ".png"
                patch_5 = subpatch.resize(size=(boxSz,boxSz), resample=Image.ANTIALIAS)
                ws_5x = self.whitespace_check(im=patch_5)
                if ws_5x < 0.9:
                    patch_5.save(os.path.join(fname, '5x', subpatch_name))

                #now we go to 10x which has 4 boxes in the 5x box
                xsubnum = 1
                for i in range(1,int(sz5x/sz10x+1)):
                    for j in range(1,int(sz5x/sz10x+1)):
                        xsubpatch = subpatch.crop(box=(sz10x*(i-1),sz10x*(j-1),sz10x*i,sz10x*j))
                        xsubpatch_name = tma_id + "_sub" + str(patchnum) + "_xsub" + str(xsubnum) + "_" + tma_outcome + ".png"
                        patch_10 = xsubpatch.resize(size=(boxSz, boxSz), resample=Image.ANTIALIAS)
                        ws_10x = self.whitespace_check(im=patch_10)
                        if ws_10x < 0.9:
                            patch_10.save(os.path.join(fname, '10x', xsubpatch_name))

                        xxsubnum = 1
                        for a in range(1,int(sz10x/boxSz+1)):
                            for b in range(1,int(sz10x/boxSz+1)):
                                xxsubpatch = xsubpatch.crop(box=(boxSz*(a-1),boxSz*(b-1),boxSz*a,boxSz*b))
                                xxsubpatch_name = tma_id + "_sub" + str(patchnum) + "_xsub" + str(xsubnum) + "_" + "_xxsub" + str(xxsubnum) + "_" + tma_outcome + ".png"
                                ws_20x = self.whitespace_check(im=xxsubpatch)
                                if ws_20x < 0.9:
                                    xxsubpatch.save(os.path.join(fname, '20x', xxsubpatch_name))
                                xxsubnum +=1

                        xsubnum += 1

                patchnum += 1

        tma_img = patch.resize(size=self.tma_size,resample=Image.ANTIALIAS)
        tma_savename = tma_id + tma_outcome + ".png"
        tma_img.save(os.path.join(fname,'tma',tma_savename))


    def tma_doSpatial(self,slide_name,startx,starty,tma_id, tma_outcome):
        oslide = openslide.OpenSlide(os.path.join(self.data_location,slide_name))
        shape = self.tma_size #width = oslide.dimensions[0]; #height = oslide.dimensions[1];
        patch = oslide.read_region((startx, starty), 0, shape);
        fname = os.path.join(self.save_location,'TMA',tma_id);
        os.mkdir(fname)
        xml_name = os.path.join(self.data_location, slide_name.replace('svs', 'xml'))
        #if xml exists:
            #tma_mask = self.read_xml(xml_file=xml_name,shape=shape)
            #make saveDir name
            #plt.imsave("C:/Users/harmonsa/Documents/mask_file.png", tma_mask.astype('uint'))


    def whitespace_check(self,im):
        bw = im.convert('L')
        bw = np.array(bw)
        bw = bw.astype('float')
        bw=bw/255
        prop_ws = (bw > 0.8).sum()/(bw>0).sum()
        return prop_ws

    def read_xml(self,xml_file,shape):
        xml = minidom.parse(xml_file)
        regions_ = xml.getElementsByTagName("Region")
        regions, region_labels = [], []
        counter = 1
        mask_all = np.zeros(shape, dtype=np.int)
        nmask_all = np.zeros(shape, dtype=np.int)
        for region in regions_:
            #print(region)
            r_label = region.getAttribute('Text')
            #print(r_label)
            type = int(region.getAttribute("NegativeROA"))
            print(str(type))
            vertices = region.getElementsByTagName("Vertex")
            xcoords = np.zeros((len(vertices), 1))
            ycoords = np.zeros((len(vertices), 1))
            coords = np.zeros((len(vertices), 2))
            for i, vertex in enumerate(vertices):
                xcoords[i][0] = vertex.attributes['X'].value
                ycoords[i][0] = vertex.attributes['Y'].value
                coords[i][0] = vertex.attributes['X'].value
                coords[i][1] = vertex.attributes['Y'].value

            regions.append(coords)
            mask = self.poly2mask(vertex_row_coords=ycoords, vertex_col_coords=xcoords, shape=shape)
            if type == 0:
                mask_all += mask
            elif type == 1:
                nmask_all += mask
            counter = counter + 1
        print('found ' + str(counter) + ' regions for parsing')
        mask_agg = np.ma.masked_where(mask_all > 0, mask_all)
        nmask_agg = np.ma.masked_where(nmask_all > 0, nmask_all)
        final_mask = (mask_agg.mask.astype(np.int) - nmask_agg.mask.astype(np.float32)).astype(np.bool)
        return final_mask


    def poly2mask(self,vertex_row_coords, vertex_col_coords, shape):
        ''''''
        fill_row_coords, fill_col_coords = draw.polygon(vertex_row_coords, vertex_col_coords, shape)
        mask = np.zeros(shape, dtype=np.int)
        mask[fill_row_coords, fill_col_coords] = 1
        return mask


    def get_binary_mask(self,shape,regions,region_labels):
        ''' creates a mask based on region coordinates, returns all coordinates contained in mask'''

        div = 1  # see sandbox note below
        maskshape = shape
        region_points = {}
        print("creating mask for " + str(len(set(region_labels))) + " label types")
        for label in set(region_labels):
            label_list = [i for i, e in enumerate(region_labels) if e == label]
            labels = tuple(regions[i] for i in label_list)
            label_points = np.empty((0,2),int)
            print("   ... parsing " + str(len(labels)) + " regions associated with this label")
            count = 1
            for reg in labels:
                print(str(count))
                xcoords = tuple(int(reg[i, 0]/div) for i in range(0, len(reg)))
                ycoords = tuple(int(reg[i, 1]/div) for i in range(0, len(reg)))
                points = self.poly2mask(vertex_row_coords=ycoords, vertex_col_coords=xcoords, shape=maskshape)
                label_points= np.append(label_points,points,axis=0)
                count+=1
            region_points[label] = label_points

        return region_points

    def poly2mask(self, vertex_row_coords, vertex_col_coords, shape):
        ''''''
        fill_row_coords, fill_col_coords = draw.polygon(vertex_row_coords, vertex_col_coords, shape)
        points = np.zeros((len(fill_row_coords), 2))
        points[:,0] = fill_col_coords
        points[:,1] = fill_row_coords
        return points


if __name__ == '__main__':
    c = TMA_processFile()
    c.organize_tmas()
    #c.read_tmaimage()
