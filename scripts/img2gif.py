import imageio
import glob
import os
import cv2
from natsort import natsorted
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--src_path', type=str, required=True)
parser.add_argument('--dst_path', type=str, required=True)
parser.add_argument('--resolution', type=int, default=200)
args = parser.parse_args()

def load_img(name, resolution):
    img = imageio.imread(filename)
    img = cv2.resize(img, (resolution, resolution), interpolation=cv2.INTER_AREA) # resize to 200x200
    return img

if __name__ == '__main__':
    filenames = glob.glob(os.path.join(args.src_path, "*.png"))
    filenames = natsorted(filenames)

    images = []
    for filename in filenames:
        images.append(load_img(filename, args.resolution))
    imageio.mimsave(args.dst_path, images)