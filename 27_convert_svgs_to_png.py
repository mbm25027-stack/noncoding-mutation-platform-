import os
import re
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


def find_browser() -> str | None:
    candidates = [
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    for name in ("msedge", "chrome"):
        found = shutil.which(name)
        if found:
            return found
    return None


def get_svg_dimensions(svg_path: Path) -> tuple[int, int]:
    text = svg_path.read_text(encoding="utf-8", errors="ignore")[:2000]
    w_match = re.search(r'\bwidth=["\'](\d+(?:\.\d+)?)', text)
    h_match = re.search(r'\bheight=["\'](\d+(?:\.\d+)?)', text)
    if w_match and h_match:
        return int(round(float(w_match.group(1)))), int(round(float(h_match.group(1))))
    vb_match = re.search(r'\bviewBox=["\'][\d.\-]+\s+[\d.\-]+\s+([\d.]+)\s+([\d.]+)', text)
    if vb_match:
        return int(round(float(vb_match.group(1)))), int(round(float(vb_match.group(2))))
    return 1080, 480


def read_png_dimensions(png_path: Path) -> tuple[int, int]:
    data = png_path.read_bytes()
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        return 0, 0
    w, h = struct.unpack(">II", data[16:24])
    return int(w), int(h)


def convert_svgs_to_png(fig_dir: Path, stems: list[str], scale: int = 3):
    browser = find_browser()
    if not browser:
        print("Error: Microsoft Edge / Chrome executable not found.")
        return

    print(f"Using system browser renderer: {browser} (scale={scale}x ~ 300 DPI)\n")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        for stem in stems:
            svg_file = fig_dir / f"{stem}.svg"
            png_file = fig_dir / f"{stem}.png"
            if not svg_file.exists():
                if png_file.exists():
                    w_px, h_px = read_png_dimensions(png_file)
                    print(f"[SKIP] {stem}.svg not found, but {stem}.png already exists ({w_px}x{h_px} px, {png_file.stat().st_size/1024:.1f} KB)")
                else:
                    print(f"[MISS] Neither {stem}.svg nor {stem}.png found in {fig_dir}")
                continue

            w, h = get_svg_dimensions(svg_file)
            html_file = tmp_path / f"{stem}.html"
            svg_uri = svg_file.resolve().as_uri()
            html_content = (
                "<!DOCTYPE html><html><head><meta charset='utf-8'>"
                "<style>*{margin:0;padding:0;box-sizing:border-box;}html,body{width:"
                f"{w}px;height:{h}px;overflow:hidden;background:#ffffff;}}"
                f"img{{display:block;width:{w}px;height:{h}px;}}</style></head>"
                f"<body><img src='{svg_uri}'/></body></html>"
            )
            html_file.write_text(html_content, encoding="utf-8")

            cmd = [
                browser,
                "--headless=new",
                "--disable-gpu",
                "--hide-scrollbars",
                f"--force-device-scale-factor={scale}",
                f"--window-size={w},{h}",
                f"--screenshot={str(png_file.resolve())}",
                html_file.resolve().as_uri(),
            ]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

            if png_file.exists() and png_file.stat().st_size > 1000:
                w_px, h_px = read_png_dimensions(png_file)
                size_kb = png_file.stat().st_size / 1024.0
                print(f"[OK]   {png_file.name:<44} -> {w_px} x {h_px} px ({size_kb:.1f} KB)")
            else:
                print(f"[FAIL] Could not render {png_file.name}")


if __name__ == "__main__":
    FIG_DIR = Path(r"D:\noncoding mutation finding surver\results\figures")
    TARGET_FIGURES = [
        "fig17_eqtl_expression_correlation",
        "fig18_ml_pr_roc_benchmark",
        "fig19_regulatory_architecture_loocv",
        "fig20_pancancer_validation_overview",
    ]
    convert_svgs_to_png(FIG_DIR, TARGET_FIGURES, scale=3)