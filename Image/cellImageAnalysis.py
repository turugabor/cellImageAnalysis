import matplotlib.pyplot as plt
from skimage import io
import numpy as np

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
            
class Experiment:
    
    def __init__(self, path):
        self.path = path
        self.images = []
        
    def get_images(self):
        pass
    
    def analyse(self):
        pass   
    