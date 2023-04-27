import matplotlib.pyplot as plt
from skimage import io
import numpy as np

from skimage.filters import threshold_mean
from skimage.morphology import binary_closing, binary_opening
from skimage.measure import label
from scipy import ndimage as ndi
from skimage.segmentation import watershed
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
    
    #az egyik órán bemutatott megoldás, amit lemásoltam
    
    def __init__(self, path, place)
    self.path = path
    self.place = place
    
    def load_ximage(self)
        tifs = [f for f in listdir(self.path) if isfile(join(self.path, f)) and f.endswith(".tif")]
        
        matches = [match for match in tifs if self.place in match]
        
        for i in range(len(matches)):
            matches[i] = join(self.path, matches[i])
            
        img = np.array(0)
        for i in matches
        
            img = img+io.imread(i)-io.imread(i).min()
        img = np.clip(img, 0, img.max()/3)
        self.image = img
        
    def display_ximage(self)
        plt.imshow(self.image, cmap='Greys')
        plt.axis("off")
        
    #a másik órán bemutatott megoldás, amit lemásoltam
    
class XpressImage2(Image)
    def __init__(self, path, pos, channel_names=None)
        super().__init__(path, channel_names)
        self.load_separated(pos)
        
    def load_separated(self, pos)
        self.selected = []
        self.files = os.listdir(self.path)
        for f in self.files:
            if pos in f:
                self.selected.append(self.path + f)
        img = io.imread_collection(self.selected)
        img = io.concatenate_images(img)
        img = np.swapaxes(img, 0, 2)
        img = np.swapaxes(img, 0, 1)
        self.image = img
            
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
    