# -*- coding: utf-8 -*-
from pathlib import Path
import re
from PIL import Image, ImageDraw, ImageFont

# ================== 你只需要改这里 ==================
CE_DIR  = Path("/mnt/data/SYC/xasfitting/newversion/zif8result")   # 以 Ce-L3 文件夹为主
LA_DIR  = Path("/mnt/data/SYC/xasfitting/full_LaCe_doublelayer/zif8/la")   # 去这里找对应 La-L2
OUT_DIR = Path("/mnt/data/SYC/xasfitting/newversion/combine") # 输出到这里
STITCH_MODE = "h"  # "h"=左右拼, "v"=上下拼
# ====================================================

IMG_EXTS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp", ".bmp"}

# 允许 Ce-L3- 或 Ce-L3_ 或 Ce-L3 直接连接；同理 La-L2
CE_RE = re.compile(r"(?:^|.*?)(?:Ce-L3[-_]?)(?P<key>.+?)_Ce_L3$", re.IGNORECASE)
LA_RE = re.compile(r"(?:^|.*?)(?:La-L2[-_]?)(?P<key>.+?)_La_L2$", re.IGNORECASE)

def is_image(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in IMG_EXTS

def key_from_ce_stem(stem: str):
    m = CE_RE.match(stem)
    return m.group("key") if m else None

def key_from_la_stem(stem: str):
    m = LA_RE.match(stem)
    return m.group("key") if m else None

def safe_font(size: int = 36) -> ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for f in candidates:
        try:
            return ImageFont.truetype(f, size=size)
        except Exception:
            pass
    return ImageFont.load_default()

def add_title_bar(img: Image.Image, title: str, bar_h: int = 80) -> Image.Image:
    w, h = img.size
    out = Image.new("RGB", (w, h + bar_h), (255, 255, 255))
    out.paste(img.convert("RGB"), (0, bar_h))

    draw = ImageDraw.Draw(out)
    font = safe_font(36)
    bbox = draw.textbbox((0, 0), title, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    x = max(10, (w - tw) // 2)
    y = (bar_h - th) // 2
    draw.text((x, y), title, fill=(0, 0, 0), font=font)
    draw.line([(0, bar_h - 1), (w, bar_h - 1)], fill=(200, 200, 200), width=2)
    return out

def stitch_two(a: Image.Image, b: Image.Image, mode: str = "h") -> Image.Image:
    a = a.convert("RGB")
    b = b.convert("RGB")

    if mode == "h":
        target_h = max(a.size[1], b.size[1])
        def resize_to_h(im):
            w, h = im.size
            if h == target_h:
                return im
            new_w = int(round(w * target_h / h))
            return im.resize((new_w, target_h), Image.LANCZOS)

        a2 = resize_to_h(a)
        b2 = resize_to_h(b)
        out = Image.new("RGB", (a2.size[0] + b2.size[0], target_h), (255, 255, 255))
        out.paste(a2, (0, 0))
        out.paste(b2, (a2.size[0], 0))
        return out

    if mode == "v":
        target_w = max(a.size[0], b.size[0])
        def resize_to_w(im):
            w, h = im.size
            if w == target_w:
                return im
            new_h = int(round(h * target_w / w))
            return im.resize((target_w, new_h), Image.LANCZOS)

        a2 = resize_to_w(a)
        b2 = resize_to_w(b)
        out = Image.new("RGB", (target_w, a2.size[1] + b2.size[1]), (255, 255, 255))
        out.paste(a2, (0, 0))
        out.paste(b2, (0, a2.size[1]))
        return out

    raise ValueError("STITCH_MODE must be 'h' or 'v'")

def build_la_index(la_dir: Path) -> dict[str, Path]:
    idx = {}
    for p in la_dir.rglob("*"):
        if not is_image(p):
            continue
        k = key_from_la_stem(p.stem)
        if not k:
            continue
        # 若同 key 多张，取路径排序最后一张（你也可改成取最新时间）
        if (k not in idx) or (str(p) > str(idx[k])):
            idx[k] = p
    return idx

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    la_index = build_la_index(LA_DIR)

    total = paired = single = skipped = 0

    ce_images = sorted([p for p in CE_DIR.rglob("*") if is_image(p)])
    for ce_img in ce_images:
        key = key_from_ce_stem(ce_img.stem)
        if not key:
            skipped += 1
            continue

        total += 1
        la_img = la_index.get(key)

        try:
            ce_im = Image.open(ce_img)
        except Exception as e:
            print(f"[WARN] cannot open Ce image: {ce_img} ({e})")
            continue

        if la_img and la_img.exists():
            try:
                la_im = Image.open(la_img)
                merged = stitch_two(ce_im, la_im, mode=STITCH_MODE)
                merged = add_title_bar(merged, key)
                out_path = OUT_DIR / f"{key}.png"
                merged.save(out_path, dpi=(300, 300))
                paired += 1
                continue
            except Exception as e:
                print(f"[WARN] failed pair {ce_img.name} + {la_img.name}: {e}")

        # 找不到 La-L2 或拼接失败：只输出 Ce-L3
        merged = add_title_bar(ce_im.convert("RGB"), key)
        out_path = OUT_DIR / f"{key}.png"
        merged.save(out_path, dpi=(300, 300))
        single += 1

    print(f"[DONE] total(valid Ce): {total}, paired: {paired}, single: {single}, skipped(non-matching): {skipped}")
    print(f"[OUT]  {OUT_DIR.resolve()}")

if __name__ == "__main__":
    main()
