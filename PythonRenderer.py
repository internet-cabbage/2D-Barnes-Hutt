import numpy as np
import vispy.scene
from vispy.scene import visuals, transforms
from tqdm import tqdm
from moviepy import VideoClip

filePath = "/Users/luthaisb/Code/C++/BarnesHutt/DataOutput.bin"

with open(filePath, 'rb') as f:
    N, tSteps = np.fromfile(f, np.int32, count=2)
    data = np.fromfile(f, np.float32).reshape(-1, N, 2)

print('Initial shape: ', data.shape)

print('N value:', N)
print('tSteps value:', tSteps)




# Make the canvas and add a simple viewer
canvas = vispy.scene.SceneCanvas(keys='interactive', show=True)
view = canvas.central_widget.add_view()

scatter = visuals.Markers(parent=view.scene)
scatter.set_data(data[0])

maxDist = np.max(data[0])
print("maxDist: ", maxDist)
view.camera = 'panzoom'
view.camera.set_range(x=(-maxDist, maxDist), y=(-maxDist, maxDist))


inc = 0
def frame(event):
    global inc
    if inc < tSteps:
        scatter.set_data(data[inc])
        inc += 1
    else:
        inc = 0
        print("Looped")
print("NaN:", np.isnan(data).any(), "  Inf:", np.isinf(data).any())
timer = vispy.app.Timer(interval='auto',connect=frame,start=True)

vispy.app.run()