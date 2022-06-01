
# Combine several GIFS -------------------------------------
# Code updated from https://stackoverflow.com/questions/51517685/combine-several-gif-horizontally-python

# Load packages
import imageio
import numpy as np 

# Imports gifs
gif1 = imageio.get_reader('bouts_anim.gif')
gif2 = imageio.get_reader('vo2_anim.gif')
gif3 = imageio.get_reader('hr_anim.gif')
gif4 = imageio.get_reader('tcpo2_anim.gif')

# Set number of frames
number_of_frames = gif1.get_length()

# Create writer object
new_gif = imageio.get_writer('output.gif')

# Create new gif
for frame_number in range(number_of_frames):
    img1 = gif1.get_next_data()
    img2 = gif2.get_next_data()
    img3 = gif3.get_next_data()
    img4 = gif4.get_next_data()
    #here is the magic
    new_image = np.vstack((img1, img2, img3, img4))
    new_gif.append_data(new_image)

gif1.close()
gif2.close()   
gif3.close()
gif4.close()
new_gif.close()

# Create mp4 file -----------------------------------------------
import os
os.environ["IMAGEIO_FFMPEG_EXE"] = "C:\\ffmpeg\\bin\\ffmpeg"

# Code updated from https://stackoverflow.com/questions/40726502/python-convert-gif-to-videomp4
import ffmpeg
import moviepy.editor as mp
clip = mp.VideoFileClip("output.gif")
clip.write_videofile("output.mp4", codec='h264_qsv', bitrate = '1100k', audio_codec='avc1', fps = 25, audio = False)
