import warnings
from tifffile import imwrite
warnings.filterwarnings("ignore")

def savetiff(file_name, volume, meta_data):
    try:
        imwrite(file_name,
                volume,
                imagej=True,
                metadata=meta_data)
    except IOError:
        status = -1
    else:
        status = 0
    
    return status

status = savetiff(file, vol, mdata)