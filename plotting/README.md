Two example scripts for reading and plotting the output in Python (assuming at least Python 3). 

Assuming you have bfield.****** and bfield.hdr in this directory (the ones I used here are too large to upload, but are a shapshot from this [paper](https://iopscience.iop.org/article/10.3847/1538-4357/ad7941)).

After executing ./fieldline in the directory, the outputs will be placed in the PFLS folder with the tracing parameters set within fieldline.cnt. The are two examples here, one for a vertical slice (plane of const. phi) and another traced from the lower boudnary (R = 1). I've renamed fieldline.cnt and PFLS with _vert and _lower for each, respectively.

Run read_qmap_vert.py or read_qmap_photo.py to plot the various quantities. The .png files show examples of what should be produced. 
