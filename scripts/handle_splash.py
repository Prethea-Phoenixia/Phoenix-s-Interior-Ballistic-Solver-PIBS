"""
cr: https://github.com/pyinstaller/pyinstaller/issues/8579
fix border of splash image being pink;
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def check_alpha(img: np.ndarray):
    contain_translucent = (img[..., 3] < 1) & (img[..., 3] > 0)
    return contain_translucent.any()


def processing_image(img_path):
    img_path = Path(img_path)
    folder = img_path.parent
    img_name = img_path.stem
    img = plt.imread(img_path)
    print("Contain translucency pixels (Before):", check_alpha(img))
    img[..., 3] = np.where(img[..., 3] != 1, 0, 1)
    print("Contain translucency pixels (After):", check_alpha(img))
    plt.imsave(folder.joinpath(f"{img_name}_out.png"), img, dpi=300)
    return


if __name__ == "__main__":
    processing_image("splash.png")
