import sys
from PIL import Image
from render import still
# usage: sheet.py name scene:t scene:t scene:t scene:t
name = sys.argv[1]; specs = sys.argv[2:]
ims = []
for sp in specs:
    i, t = sp.split(":"); ims.append(Image.open(still(int(i) - 1, float(t))))
sheet = Image.new("RGB", (1920, 1080 * ((len(ims) + 1) // 2) // 2))
for k, im in enumerate(ims):
    sheet.paste(im.resize((960, 540), Image.LANCZOS), ((k % 2) * 960, (k // 2) * 540))
sheet.save(f"out/stills/{name}.png"); print(f"out/stills/{name}.png")
