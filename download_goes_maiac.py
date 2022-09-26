import re, subprocess
from datetime import date, timedelta
from pathlib import Path
import requests

target_dir = Path('/mnt/qnap_geo/goes16_maiac')
base_url = 'https://data.nas.nasa.gov/geonex/geonexdata/GOES16/GEONEX-L2/MAIAC'

def seq(a, b):
    return range(a, b + 1)

def GET(*args, **kwargs):
    r = requests.get(*args, **kwargs)
    r.raise_for_status()
    return r.text

conus_tiles = sorted(
  # http://web.archive.org/web/20220923195042/https://www.nasa.gov/sites/default/files/thumbnails/image/globalgrid_v3.png
    f'h{h:02}v{v:02}'
    for v, hs in {
        1: seq(9, 18),
        2: seq(9, 18),
        3: seq(9, 18),
        4: seq(9, 17),
        5: seq(11, 16)}.items()
    for h in hs
    if not (
       # Not currently present on the server.
       v == 1 or
       h == 18))

start, end = date(2018, 1, 1), date(2019, 12, 31)
tile_days = [(tile, start + timedelta(days = i))
    for i in seq(0, (end - start).days)
    for tile in conus_tiles]

for i, (tile, day) in enumerate(tile_days):
    print(f'{tile} {day} - {i + 1 :,} / {len(tile_days):,}')
    out_dir = target_dir / str(day.year) / tile
    out_dir.mkdir(parents = True, exist_ok = True)
    dir_url = f'{base_url}/{tile}/{day.year}/{day.strftime("%j")}/'
    for file_url in re.findall('<td><a href="([^"]+)', GET(dir_url)):
        file_url = requests.compat.urljoin(dir_url, file_url)
        subprocess.run(('curl',
            '--fail', '--remote-time',
            str(file_url),
            '-o', str(out_dir / Path(file_url).name)))
