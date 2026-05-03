from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"C:\TRKB_WP2\docking")
VIEW_DIR = ROOT / "figure_batch" / "views"
PANEL_DIR = ROOT / "figure_batch" / "panels"
PANEL_DIR.mkdir(parents=True, exist_ok=True)

VIEW_NAMES = ["front", "front45", "side", "back"]

def load_font(size):
    candidates = [
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibri.ttf",
        r"C:\Windows\Fonts\arialbd.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size)
    return ImageFont.load_default()

title_font = load_font(52)
label_font = load_font(34)

# 自动找所有 label
all_pngs = list(VIEW_DIR.glob("*_front.png"))
labels = [p.name.replace("_front.png", "") for p in all_pngs]

for label in labels:
    images = []
    for vn in VIEW_NAMES:
        img_path = VIEW_DIR / f"{label}_{vn}.png"
        if not img_path.exists():
            print(f"Missing {img_path}")
            break
        images.append(Image.open(img_path).convert("RGB"))
    else:
        # 统一尺寸
        w = min(im.width for im in images)
        h = min(im.height for im in images)
        images = [im.resize((w, h)) for im in images]

        margin = 40
        header_h = 120
        caption_h = 60
        panel_w = w * 2 + margin * 3
        panel_h = header_h + (h + caption_h) * 2 + margin * 3

        canvas = Image.new("RGB", (panel_w, panel_h), "white")
        draw = ImageDraw.Draw(canvas)

        # 标题
        title = label
        tw = draw.textbbox((0, 0), title, font=title_font)[2]
        draw.text(((panel_w - tw) // 2, 25), title, fill="black", font=title_font)

        positions = [
            (margin, header_h + margin),
            (w + margin * 2, header_h + margin),
            (margin, header_h + h + caption_h + margin * 2),
            (w + margin * 2, header_h + h + caption_h + margin * 2),
        ]

        captions = ["Front", "45°", "Side", "Back"]

        for img, pos, cap in zip(images, positions, captions):
            x, y = pos
            canvas.paste(img, (x, y))
            bbox = draw.textbbox((0, 0), cap, font=label_font)
            cap_w = bbox[2] - bbox[0]
            draw.text((x + (w - cap_w) // 2, y + h + 10), cap, fill="black", font=label_font)

        out_path = PANEL_DIR / f"{label}_panel.png"
        canvas.save(out_path, dpi=(300, 300))
        print(f"Saved {out_path}")

print("All panels done.")