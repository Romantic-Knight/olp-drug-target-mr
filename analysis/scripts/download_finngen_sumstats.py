# Download FinnGen R13 OLP sumstats (full files) with resume support
import urllib.request, os, time, sys

BASE = 'https://storage.googleapis.com/finngen-public-data-r13/summary_stats/'
DEST_DIR = r'D:\AI-Wisp-Science\Drug Target MR\data\raw\finngen_r13'
FILES = [
    'finngen_R13_K11_LICHEN_PLANUS_WIDE.gz',
    'finngen_R13_K11_ORAL_LICHEN_PLANUS_INCLAVO.gz',
]

def get_size(url):
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get('Content-Length', 0))

def download(url, dest, expected):
    """Resume-capable download. Returns (ok, size)."""
    if os.path.exists(dest) and os.path.getsize(dest) == expected:
        print(f'  already complete: {dest}', flush=True)
        return True, expected
    partial = dest + '.part'
    have = os.path.getsize(partial) if os.path.exists(partial) else 0
    headers = {'User-Agent': 'Mozilla/5.0'}
    if have:
        headers['Range'] = f'bytes={have}-'
    req = urllib.request.Request(url, headers=headers)
    t0 = time.time()
    last = have
    last_t = t0
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            mode = 'ab' if have else 'wb'
            with open(partial, mode) as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    now = time.time()
                    if now - last_t >= 10:
                        spd = (os.path.getsize(partial) - last) / (now - last_t) / 1e6
                        print(f'  {now-t0:.0f}s: {os.path.getsize(partial)/1e6:.0f}/{expected/1e6:.0f} MB  ({spd:.1f} MB/s)', flush=True)
                        last = os.path.getsize(partial); last_t = now
        sz = os.path.getsize(partial)
        if sz == expected:
            os.replace(partial, dest)
            print(f'  DONE {os.path.basename(dest)}: {sz/1e6:.1f} MB in {time.time()-t0:.0f}s', flush=True)
            return True, sz
        else:
            print(f'  partial {sz/1e6:.1f}/{expected/1e6:.1f} MB (connection closed; will resume next run)', flush=True)
            return False, sz
    except Exception as e:
        print(f'  error: {e} (partial kept for resume)', flush=True)
        return False, os.path.getsize(partial) if os.path.exists(partial) else 0

os.makedirs(DEST_DIR, exist_ok=True)
for f in FILES:
    url = BASE + f
    dest = os.path.join(DEST_DIR, f)
    expected = get_size(url)
    print(f'== {f}: {expected/1e6:.1f} MB', flush=True)
    ok, sz = download(url, dest, expected)
    if not ok:
        print(f'  WARNING: not complete ({sz/1e6:.1f}/{expected/1e6:.1f} MB)', flush=True)
print('all done', flush=True)
