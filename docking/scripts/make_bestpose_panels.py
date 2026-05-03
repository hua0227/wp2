from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"C:\TRKB_WP2\docking")
VIEW_DIR = ROOT / "figure_best" / "views"
PANEL_DIR = ROOT / "figure_best" / "panels"
PANEL_DIR.mkdir(parents=True, exist_ok=True)

VIEW_NAMES = ["front", "front45", "side", "back"]

def load_font(size):
    candidates = [
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\arialbd.ttf",
        r"C:\Windows\Fonts\calibri.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size)
    return ImageFont.load_default()

title_font = load_font(46)
label_font = load_font(30)

front_files = sorted(VIEW_DIR.glob("*_front.png"))
labels = [p.name.replace("_front.png", "") for p in front_files]

print("Panels to make:", len(labels))

for label in labels:
    imgs = []
    ok = True
    for vn in VIEW_NAMES:
        img_path = VIEW_DIR / f"{label}_{vn}.png"
        if not img_path.exists():
            print("Missing:", img_path)
            ok = False
            break
        imgs.append(Image.open(img_path).convert("RGB"))
    if not ok:
        continue

    w = min(im.width for im in imgs)
    h = min(im.height for im in imgs)
    imgs = [im.resize((w, h)) for im in imgs]

    margin = 35
    header_h = 100
    caption_h = 50
    panel_w = w * 2 + margin * 3
    panel_h = header_h + (h + caption_h) * 2 + margin * 3

    canvas = Image.new("RGB", (panel_w, panel_h), "white")
    draw = ImageDraw.Draw(canvas)

    title = label
    bbox = draw.textbbox((0, 0), title, font=title_font)
    tw = bbox[2] - bbox[0]
    draw.text(((panel_w - tw) // 2, 20), title, fill="black", font=title_font)

    positions = [
        (margin, header_h + margin),
        (w + margin * 2, header_h + margin),
        (margin, header_h + h + caption_h + margin * 2),
        (w + margin * 2, header_h + h + caption_h + margin * 2),
    ]
    captions = ["Front", "45掳", "Side", "Back"]

    for img, pos, cap in zip(imgs, positions, captions):
        x, y = pos
        canvas.paste(img, (x, y))
        bb = draw.textbbox((0, 0), cap, font=label_font)
        cw = bb[2] - bb[0]
        draw.text((x + (w - cw) // 2, y + h + 8), cap, fill="black", font=label_font)

    out_path = PANEL_DIR / f"{label}_panel.png"
    canvas.save(out_path, dpi=(300, 300))
    print("Saved:", out_path)

print("All panels done.")
