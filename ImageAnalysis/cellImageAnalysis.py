import matplotlib.pyplot as plt
from skimage import io
import numpy as np
import pandas as pd

from skimage.filters import threshold_mean
from skimage.morphology import binary_closing, binary_opening
from skimage.measure import label
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from os import listdir
from os.path import isfile, join

from skimage.morphology import binary_erosion

from cellpose import models

class Image:
    
    def __init__(self, path, channel_names=None):
        self.path = path
        self.channel_names = channel_names
        self.name = ""
        
    def load_image(self):
        img = io.imread(self.path)
        
        if img.shape[0] <6:
            img = np.swapaxes(img, 0, 2)
            img = np.swapaxes(img, 0, 1)
            
        if img.shape[1] <6:
            img = np.swapaxes(img, 1, 2)
        
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
    
#     def save(self, path):
#         io.imsave(path+"/"+self.place+".tif", self.image)
    
class XpressImage(Image):

    
    def __init__(self, path, place, channel_names=None):
        super().__init__(path)
        self.place = place
    
    def load_image(self):   #betölt és összerak
        #leszűrjük, hogy a mappából csak a tif file-okat vizsgáljuk
        tifs = [f for f in listdir(self.path) if isfile(join(self.path, f)) and  f.endswith(".tif")]
        
        #majd kiválasztjuk belőle a megadott hely alapján a megfelelő fileokat
        matches = [match for match in tifs if self.place in match]
        
        #teljes elérési utat hozunk létre a fileokhoz, 
        for i in range(len(matches)):
            matches[i]=join(self.path,matches[i])
        
        
        #img megkapja az elérési utakat, összerakja az array-t
        #/65535*255
        img = io.imread_collection(matches)
        img = (np.stack(img, axis=2))
        
               
        self.image=img
               
    def display_ximage(self):
        """
        Az Ximage esetén létrehoz a 3 csatornából egy képet. Elsősorban debugging céljából került be.
        """
        #létrehozunk egy arrayt a kép dimenzióival
        dim1=self.image.shape[0]
        dim2=self.image.shape[1]
        img=np.zeros((dim1,dim2))
        
        #a csatornákat összeadjuk
        channels = self.image.shape[-1]
        for i in range(channels):
            #mindenhonnan levonjuk a hátteret egyesével, mivel a külön csatornákon eltérhet
            img += self.image[:,:,i] - self.image[:,:,i].min() 
            #img=np.clip(img, 0, img.max()/3)
                
        plt.imshow(img, cmap="Greys")
        plt.axis("off")  

            
class Detector:
    
    def __init__(self, cell_channel=None, nucleus_channel=None):
        self.cell_channel = cell_channel
        self.nucleus_channel = nucleus_channel
               
    def detect_cells(self, image):
        
        image_data = image[:,:,self.cell_channel]
        thr = threshold_mean(image_data)
        binary = image_data > thr
        binary = binary_closing(binary, np.ones(shape=(8,8)))
        binary = binary_opening(binary, np.ones(shape=(15,15)))
        distance = ndi.distance_transform_edt(binary)
        coords = peak_local_max(distance, footprint=np.ones((150, 150)), labels=binary)
        mask = np.zeros(distance.shape, dtype=bool)
        mask[tuple(coords.T)] = True
        markers, _ = ndi.label(mask)
        labels = watershed(-distance, markers, mask=binary)
        
        return labels
    
    def detect_nuclei(self, image):
        pass
    
    def detect_cell_borders(self, masks):
        """
        Args: 
            masks: masked image
        Returns:
            countour of the cells in the masks image"""
        
        eroded = binary_erosion(masks) 
        return (masks>0)^eroded
    
class CellposeDetector(Detector):
    
    def __init__(self, cell_channel=None, nucleus_channel=None):
        """Detects cell using the cellpose library
        Args:
            cell_channel (int): index of the cell channel in the images using zero indexing
            nucleus_channel (int): index of the nuclear channel in the images using zero indexing"""
        super().__init__(cell_channel, nucleus_channel)
        
       
        self.channels = [0,0]
        
        if self.cell_channel != None:
            self.channels[0] = self.cell_channel + 1
        if self.nucleus_channel != None:
            self.channels[1] = self.nucleus_channel + 1
      
    #gpu TRUE!!!      
        self.cell_model = models.Cellpose(gpu=True, model_type='cyto')
        self.nucleus_model = models.Cellpose(gpu=True, model_type='nuclei')

        
    def detect_cells(self, image, diameter=100):
        cell_masks, __, __, __ = self.cell_model.eval(image, diameter=diameter, channels=self.channels)
        return cell_masks
    
    def detect_nuclei(self, image, diameter=100):
        nuclear_masks, __, __, __ = self.nucleus_model.eval(image, diameter=diameter, channels=self.channels)
        return nuclear_masks

    

class Analyzer:
    
    def __init__(self):
        pass
    
    def get_cell_fluorescence(self, img, masks):
        img = img.reshape(-1,1)
        masks = masks.reshape(-1,1)
        
        joint = np.concatenate((img, masks), axis=1)
        
        result = pd.DataFrame(joint)
        
        result.columns = ["intensity", "mask"]
        
        return result.groupby("mask").mean()

    
class Experiment:
    """
    Uses other classes (written for this demo) to tell the fluorescence intensity of cells.
    Needs only the path to the folder with pictures.
    
    """
    
    def __init__(self, path, cell_channel=None, nucleus_channel=None, signal_channel=None):
        self.path = path
        self.cell_channel = cell_channel 
        self.signal_channel = signal_channel
        self.nucleus_channel = nucleus_channel
        
        #meghívjuk a detectort és az analizert
        self.detector = CellposeDetector(self.cell_channel, self.nucleus_channel)
        self.analyzer = Analyzer()
        
        #dataframe az eredmények gyűjtésére
        self.result = pd.DataFrame()
                        
        
    def get_images(self):
        #összeszedjük a képeket a mappából; jól bevált list comprehension
        self.images = [f for f in listdir(self.path) if isfile(join(self.path, f)) and  f.endswith(".tif")]
        print (len(self.images), "db képet találtunk, ezek: ", self.images)
    
    
    
    def analyse (self):
        #végigmegy a képeken
        for i in self.images:
            name = join(self.path, i )
                        
            #mivel a ciklusban van számítás elrejtve sok, jelezgetjük, hogy történik valami
            print(name, "loaded")
            
            #betölt
            a = Image(name)
            a.load_image()
            #csatornákat összerakjuk, így határozunk meg intenzitást
#             cell = a.image[:,:,0] + a.image[:,:,1]
            
            #itt kezdődik az "analízis";
            #sejtmaszkok
            masks = self.detector.detect_cells(a.image)
            
            #intenzitás értékek kiszámolása
            #berakjuk az eredményeket a tárolóba
            self.result[i] = self.analyzer.get_cell_fluorescence(a.image[:,:,self.signal_channel], masks)
            
            #mivel a ciklusban van számítás elrejtve sok, jelezgetjük, hogy történik valami
            print(name, "processed")
        
        print ("Process ended! You can get the results with '.result'.")
        
    def plot_hist(self):
        pass
    
   
    