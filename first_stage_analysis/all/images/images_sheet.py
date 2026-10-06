"""軌道上画像の一覧（サムネイル表）を作る．

1次データ（data/downlink/<sat>/images.csv とそこに列挙された画像）だけを入力とし，
加工は縮小と並べ替えのみ．

入力
  data/downlink/<sat>/images.csv   画像の一覧（satellite, instrument, mission, path, status, …）
  data/downlink/<sat>/images/<instrument>/*   画像本体
出力
  first_stage_analysis/<sat>/images/images_<instrument>[_<mission>].jpg
    衛星×カメラごとに 1 枚．複数ミッションを持つカメラはミッションごとに分ける．
    撮影日時（ファイル名）順に並べる．部分受信で末尾が欠けた画像は復号できた範囲まで表示し，
    復号できない画像は灰色の枠に images.csv の status を添えて示す．
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFile, ImageFont

ImageFile.LOAD_TRUNCATED_IMAGES = True  # 部分受信の画像を復号できた範囲まで読む

ROOT = Path(__file__).resolve().parents[3]
DL = ROOT / "data" / "downlink"
FIRST = ROOT / "first_stage_analysis"  # 出力は first_stage_analysis/<sat>/images/

COLS = 4
THUMB = (200, 150)      # サムネイルの最大寸法
PAD = 12                # セル間の余白
LABEL_H = 16            # ファイル名の行
TITLE_H = 28
BG = "white"
INK = "#0b0b0b"
INK2 = "#52514e"
PLACEHOLDER = "#9a9a96"


def font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.load_default(size=size)
    except TypeError:  # Pillow < 10.1
        return ImageFont.load_default()


def thumbnail(path: Path) -> Image.Image | None:
    try:
        im = Image.open(path)
        im.load()
        im = im.convert("RGB")
    except Exception:
        return None
    im.thumbnail(THUMB)
    return im


def make_sheet(title: str, rows: pd.DataFrame, out: Path) -> None:
    rows = rows.sort_values("path").reset_index(drop=True)
    n = len(rows)
    nrow = (n + COLS - 1) // COLS
    cell_w, cell_h = THUMB[0] + PAD, THUMB[1] + LABEL_H + PAD
    W = PAD + COLS * cell_w
    H = TITLE_H + nrow * cell_h + PAD
    sheet = Image.new("RGB", (W, H), BG)
    draw = ImageDraw.Draw(sheet)
    f_title, f_label = font(13), font(10)
    draw.text((PAD, 8), f"{title}  n={n}", fill=INK, font=f_title)
    for i, r in rows.iterrows():
        cx = PAD + (i % COLS) * cell_w
        cy = TITLE_H + (i // COLS) * cell_h
        im = thumbnail(ROOT / r.path)
        if im is None:
            draw.rectangle([cx, cy, cx + THUMB[0] - 1, cy + THUMB[1] - 1], fill=PLACEHOLDER)
            draw.text((cx + 6, cy + THUMB[1] // 2 - 6), f"not decodable (status: {r.status})", fill="white", font=f_label)
        else:
            sheet.paste(im, (cx + (THUMB[0] - im.width) // 2, cy + (THUMB[1] - im.height) // 2))
        # 接頭辞（衛星・ミッション）はタイトルと重複するので，日時とハッシュだけを示す
        label = "_".join(Path(r.path).stem.split("_")[-2:])
        draw.text((cx, cy + THUMB[1] + 2), label, fill=INK2, font=f_label)
    sheet.save(out, quality=85)
    print("wrote", out.relative_to(ROOT), n)


def main():
    for csv in sorted(DL.glob("*/images.csv")):
        d = pd.read_csv(csv, dtype=str).fillna("")
        for (sat, inst), g in d.groupby(["satellite", "instrument"], sort=True):
            missions = sorted(m for m in g.mission.unique())
            split = len(missions) > 1
            for m in (missions if split else [""]):
                rows = g[g.mission == m] if split else g
                title = " ".join(x for x in (sat, inst, m) if x)
                name = "_".join(x.lower() for x in ("images", inst, m) if x)
                out_dir = FIRST / csv.parent.name / "images"
                out_dir.mkdir(parents=True, exist_ok=True)
                make_sheet(title, rows, out_dir / f"{name}.jpg")


if __name__ == "__main__":
    main()
