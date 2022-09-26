import multiprocessing
import re, subprocess
import datetime as DT
from pathlib import Path
import requests

n_workers = 15
target_dir = Path('/mnt/qnap_geo/goes16_maiac')
log_path = Path('/scratch/goes16_maiac_gotten.csv')
base_url = 'https://data.nas.nasa.gov/geonex/geonexdata/GOES16/GEONEX-L2/MAIAC'

def seq(a, b):
    return range(a, b + 1)

def GET(*args, **kwargs):
    r = requests.get(*args, **kwargs)
    r.raise_for_status()
    return r.text

def get_tile_day(tile_day):
    tile, date = tile_day
    out_dir = target_dir / str(date.year) / tile
    out_dir.mkdir(parents = True, exist_ok = True)
    dir_url = f'{base_url}/{tile}/{date.year}/{date.strftime("%j")}/'
    for file_url in re.findall('<td><a href="([^"]+)', GET(dir_url)):
        file_url = requests.compat.urljoin(dir_url, file_url)
        subprocess.run(('curl',
            '--fail', '--remote-time', '--silent',
            str(file_url),
            '-o', str(out_dir / Path(file_url).name)))
    return tile_day

def main():

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

    start, end = DT.date(2018, 1, 1), DT.date(2019, 12, 31)
    tile_days = [(tile, date)
        for i in seq(0, (end - start).days)
        for date in [start + DT.timedelta(days = i)]
        for tile in conus_tiles]
    n_total = len(tile_days)
    gotten = {tuple(x.split(','))
        for x in log_path.read_text().splitlines()}
    n_gotten = len(gotten)
    tile_days = [(tile, date) for
        tile, date in tile_days
        if (str(date), tile) not in gotten]

    def status():
        print('Tile-days remaining:', format(n_total - n_gotten, ","))

    with multiprocessing.Pool(n_workers) as pool, log_path.open('at') as o_log:
        status()
        for tile, date in pool.imap_unordered(get_tile_day, tile_days):
            print(f'{date},{tile}', file = o_log)
            o_log.flush()
            n_gotten += 1
            status()

if __name__ == '__main__':
    main()
