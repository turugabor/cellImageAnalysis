import matplotlib.pyplot as plt
from skimage import io
import numpy as np

from skimage.filters import threshold_mean
from skimage.morphology import binary_closing, binary_opening
from skimage.measure import label
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from os import listdir
from os.path import isfile, join

class Image:
    
    def __init__(self, path, channel_names=None):
        self.path = path
        self.channel_names = channel_names
        
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
            plt.imshow(self.image[:,:,i], cmap="Greys") 
            if self.channel_names != None:
                plt.title(self.channel_names[i])
            plt.axis("off")
            
    def remove_background(self):
        pass
    
    
class XpressImage(Image):

    
    def __init__(self, path, place, channel_names=None):
        super().__init__(path)
        self.place = place
    
    def load_ximage(self):   #betölt és összerak
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
        
        #jó ötletnek tűntek a következőek, de ha nem teszem bele, akkor simán működik minden.
        #levonjuk a minimum értékeket, leklippeljük 0-255-re, és konvertálunk integerré (kerekítve)
        #mindezt külön-külön a csatornákon, mert nagy a szórás a min-max értékeikben        
        
        #channels=img.shape[-1]

        #for i in range(channels):
         #   img[:,:,i]= np.rint ( (np.clip( (img[:,:,i] - img[:,:,i].min() ),0,255) ) ).astype(int)

        
        self.image=img
        
    
    def display_image(self):
        super().display_image()
        
    def display_ximage(self):

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
    
    
    def save(self, path):
        io.imsave(path+"/"+self.place+".tif", self.image)

            
class Detector:
    
    def __init__(self, cell_channel=None, nucleus_channel=None):
        self.cell_channel = cell_channel
        self.nucleus_channel = nucleus_channel
        
        
    def watershed(self, binary):
        pass
    
    def make_binary(self, image):
        pass
        
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
    
    def detect_cell_borders(self, image):
        pass            
            
class Experiment:
    
    def __init__(self, path):
        self.path = path
        self.images = []
        
    def get_images(self):
        pass
    
    def analyse(self):
        pass   
    