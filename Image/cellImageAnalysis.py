import os
import matplotlib.pyplot as plt
from skimage import io
import numpy as np

from skimage.filters import threshold_mean, threshold_isodata
from skimage.morphology import binary_closing, binary_opening
from skimage.measure import label
from scipy import ndimage as ndi
from skimage.segmentation import watershed, find_boundaries
from skimage.feature import peak_local_max

class Image:
    
    def __init__(self, path, channel_names=None):
        self.path = path
        self.channel_names = channel_names
        
    def load_image(self):
        img = io.imread(self.path)
        img = np.swapaxes(img, 0, 2)
        img = np.swapaxes(img, 0, 1)
        self.image = img
        
    def display_image(self):
        channels = self.image.shape[-1]
        for i in range(channels):
            plt.subplot(1, channels, i+1)
            plt.imshow(self.image[:,:,i], cmap="gray") 
            if self.channel_names != None:
                plt.title(self.channel_names[i])
            plt.axis("off")
            
    def remove_background(self):
        pass
    
class XpressImage(Image):
    
    def __init__(self, path, name):
        if path.endswith('/')==False:
            path=path+'/'
        self.path = path
        self.imfiles=[]
        self.channel_names=[]
        for file in os.listdir(path):
            if file.startswith(name):
                self.imfiles.append(file)
                self.channel_names.append(file[len(name)+1:len(name)+3])
                

    def load_image(self):
        ims=[]
        for file in self.imfiles:
            ims.append(io.imread(self.path+file))
        self.image=np.stack(ims, axis=2)
            
class Detector:
    
    def __init__(self, cell_channel=None, nucleus_channel=None):
        self.cell_channel = cell_channel
        self.nucleus_channel = nucleus_channel
        
        
    def watershed(self, binary):
        pass
    
    def make_binary(self, image):
        pass
        
    def detect_cells(self, image):
        
        image_data = image.image[:,:,self.cell_channel]
        thr = threshold_mean(image_data)
        binary = image_data > thr
        binary = binary_closing(binary, np.ones(shape=(8,8)))
        binary = binary_opening(binary, np.ones(shape=(15,15)))
        distance = ndi.distance_transform_edt(binary)
        coords = peak_local_max(distance, footprint=np.ones((150, 150)), labels=binary)
        mask = np.zeros(distance.shape, dtype=bool)
        mask[tuple(coords.T)] = True
        markers, cnum = ndi.label(mask)
        labels = watershed(-distance, markers, mask=binary)
        
        image.cell_number=cnum
        image.cell_labels=labels
        #return labels
    
    def detect_nuclei(self, image):
        try:
            labels=image.cell_labels
        except AttributeError:
            #print('Need cell detection first.')
            #return
            self.detect_cells(image)
            labels=image.cell_labels
        cell_num=image.cell_number        
        nuclei = np.zeros(labels.shape, dtype=int)

        for cnum in range(1,cell_num,1):
            difim=image.image[:,:,self.nucleus_channel]-image.image[:,:,self.cell_channel]
            #difim=image.image[:,:,1]-image.image[:,:,0]
            difim[np.where(labels != cnum)]=0
            thr = threshold_isodata(difim)
            binary = difim > thr
            
            binary = binary_opening(binary, np.ones(shape=(4,4)))
            binary = binary_closing(binary, np.ones(shape=(4,4)))
            #binary = binary_opening(binary, np.ones(shape=(5,5)))
            binary = binary_opening(binary, np.ones(shape=(16,16)))

            nuclei[np.where(binary==True)]=cnum
            print('.', end='', flush=True) 
            
        image.nucleus_labels=nuclei
    
    def detect_cell_borders(self, image):
        try:
            labels=image.cell_labels
        except AttributeError:
            #print('Need cell detection first.')
            #return
            self.detect_cells(image)
            labels=image.cell_labels
        
        borders = find_boundaries(labels)
        image.cell_borders=borders
        #return borders
            
class Experiment:
    
    def __init__(self, path):
        self.path = path
        self.images = []
        
    def get_images(self):
        pass
    
    def analyse(self):
        pass   
    