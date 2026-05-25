#!/usr/bin/env python3
"""Generate resources/icon.ico — DNA chromatogram themed app icon.

Run once from the project root:
    python make_icon.py
"""
from __future__ import annotations
import math, struct, io
from pathlib import Path
from PIL import Image, ImageDraw


# Dark-mode GATC palette
_COLORS = {
    "G": (221, 215, 10),
    "A": (51,  210, 80),
    "T": (255, 75,  75),
    "C": (80,  148, 255),
}
_BG = (14, 20, 38)

_PEAKS = [
    ("G", 0.28, 0.88),
    ("A", 0.44, 0.72),
    ("T", 0.60, 0.95),
    ("C", 0.76, 0.65),
]
_SIGMA = 0.095


def _gauss(x: float, mu: float, sigma: float) -> float:
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _draw_icon(size: int) -> Image.Image:
    result = Image.new("RGBA", (size, size), (0, 0, 0, 0))

    pad = max(1, size // 18)

    # Background circle
    bg = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    ImageDraw.Draw(bg).ellipse(
        [pad, pad, size - pad - 1, size - pad - 1],
        fill=(*_BG, 255),
    )
    result = Image.alpha_composite(result, bg)

    y_floor = int(size * 0.84)
    max_h   = int(size * 0.64)

    for base, mu, rel_h in _PEAKS:
        color = _COLORS[base]
        layer = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        draw  = ImageDraw.Draw(layer)

        x_lo = max(0, int((mu - 3.2 * _SIGMA) * size))
        x_hi = min(size - 1, int((mu + 3.2 * _SIGMA) * size))

        poly: list[tuple[int, int]] = [(x_lo, y_floor)]
        for px in range(x_lo, x_hi + 1):
            h = _gauss(px / size, mu, _SIGMA) * max_h * rel_h
            poly.append((px, y_floor - int(h)))
        poly.append((x_hi, y_floor))

        if len(poly) >= 3:
            draw.polygon(poly, fill=(*color, 210))

        result = Image.alpha_composite(result, layer)

    # Baseline rule
    if size >= 24:
        bl = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        ImageDraw.Draw(bl).line(
            [(int(size * 0.10), y_floor + 1),
             (int(size * 0.90), y_floor + 1)],
            fill=(180, 180, 180, 140),
            width=max(1, size // 48),
        )
        result = Image.alpha_composite(result, bl)

    # Clip to circle
    mask = Image.new("L", (size, size), 0)
    ImageDraw.Draw(mask).ellipse(
        [pad, pad, size - pad - 1, size - pad - 1],
        fill=255,
    )
    result.putalpha(mask)
    return result


def _build_ico(images: list[Image.Image]) -> bytes:
    """Manually build an ICO file from a list of RGBA images.

    Pillow's ICO saver resizes a single source image; this function
    embeds each pre-rendered size as its own PNG entry.
    """
    n = len(images)
    # ICO header: 6 bytes
    # Directory entries: n * 16 bytes
    dir_offset = 6 + n * 16
    png_bufs: list[bytes] = []
    for img in images:
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        png_bufs.append(buf.getvalue())

    out = io.BytesIO()
    # ICONDIR header
    out.write(struct.pack("<HHH", 0, 1, n))  # reserved=0, type=1 (ICO), count=n

    offset = dir_offset
    for img, png in zip(images, png_bufs):
        w, h = img.size
        w8 = w if w < 256 else 0   # 0 means 256 in ICO spec
        h8 = h if h < 256 else 0
        out.write(struct.pack(
            "<BBBBHHII",
            w8, h8,      # width, height
            0, 0,        # color count, reserved
            1,           # color planes
            32,          # bits per pixel
            len(png),    # size of image data
            offset,      # offset in file
        ))
        offset += len(png)

    for png in png_bufs:
        out.write(png)

    return out.getvalue()


def main() -> None:
    out = Path("resources/icon.ico")
    out.parent.mkdir(exist_ok=True)

    sizes  = [16, 24, 32, 48, 64, 128, 256]
    images = [_draw_icon(s) for s in sizes]

    data = _build_ico(images)
    out.write_bytes(data)
    print(f"Saved {out}  ({len(data):,} bytes,  sizes={sizes})")


if __name__ == "__main__":
    main()
