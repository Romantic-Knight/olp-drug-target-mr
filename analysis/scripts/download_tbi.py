# Download FinnGen R13 .tbi index files via urllib (curl SSL fails on this host)
import urllib.request, os, time

base = 'https://storage.googleapis.com/finngen-public-data-r13/summary_stats/'
dest_dir = r'D:\AI-Wisp-Science\Drug Target MR\data\raw\finngen_r13'
os.makedirs(dest_dir, exist_ok=True)
files = ['finngen_R13_K11_LICHEN_PLANUS_WIDE.gz.tbi',
         'finngen_R13_K11_ORAL_LICHEN_PLANUS_INCLAVO.gz.tbi']

for f in files:
    url = base + f
    dest = os.path.join(dest_dir, f)
    t0 = time.time()
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest, 'wb') as out:
                total = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    total += len(chunk)
            print(f'OK {f}: {total/1e6:.2f} MB in {time.time()-t0:.1f}s', flush=True)
            break
        except Exception as e:
            print(f'  attempt {attempt+1} failed: {e}', flush=True)
            time.sleep(5)
    else:
        raise SystemExit(f'failed to download {f}')
