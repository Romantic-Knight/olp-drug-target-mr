# Zero-dependency SigV4 S3 client for deCODE summary-data downloads.
# Endpoint: https://s3-ext.decode.is (HTTPS, path-style), region us-east-1.
# Usage:
#   python s3_decode.py list [--prefix P]
#   python s3_decode.py get <key> <dest_file>
#   python s3_decode.py sync <dest_dir> [--prefix P]
import hashlib, hmac, os, sys, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from datetime import datetime, timezone

HOST = "s3-ext.decode.is"
REGION = "us-east-1"
SERVICE = "s3"
BUCKET = "largescaleplasma-2023"
ACCESS_KEY = "V1IQRNX6GFAH9RL9CF2J"
SECRET_KEY = "17LeuimS93+l9bBkmrFftiZC9quEJ1HtzEf7/Bu2"
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()


def sign(key, msg):
    return hmac.new(key, msg.encode(), hashlib.sha256).digest()


def amz_date_now():
    now = datetime.now(timezone.utc)
    return now.strftime("%Y%m%dT%H%M%SZ"), now.strftime("%Y%m%d")


def canonical_uri_for(path):
    # path = "/bucket/key..." ; quote each segment, keep slashes
    return urllib.parse.quote(path, safe="/~")


def build_auth(method, path, query_params, headers):
    amz_date, date = amz_date_now()
    payload_hash = EMPTY_SHA256
    all_headers = {
        "host": HOST,
        "x-amz-content-sha256": payload_hash,
        "x-amz-date": amz_date,
        **headers,
    }
    canonical_headers = "".join(f"{k.lower()}:{str(v).strip()}\n" for k, v in sorted(all_headers.items()))
    signed_headers = ";".join(sorted(k.lower() for k in all_headers))
    canonical_query = "&".join(
        f"{urllib.parse.quote(k, safe='-_.~')}={urllib.parse.quote(str(v), safe='-_.~')}"
        for k, v in sorted(query_params.items())
    )
    canonical_request = "\n".join([method, canonical_uri_for(path), canonical_query,
                                   canonical_headers, signed_headers, payload_hash])
    scope = f"{date}/{REGION}/{SERVICE}/aws4_request"
    string_to_sign = "\n".join(["AWS4-HMAC-SHA256", amz_date, scope,
                                hashlib.sha256(canonical_request.encode()).hexdigest()])
    k_date = sign(("AWS4" + SECRET_KEY).encode(), date)
    k_region = sign(k_date, REGION)
    k_service = sign(k_region, SERVICE)
    k_signing = sign(k_service, "aws4_request")
    signature = hmac.new(k_signing, string_to_sign.encode(), hashlib.sha256).hexdigest()
    auth = (f"AWS4-HMAC-SHA256 Credential={ACCESS_KEY}/{scope}, "
            f"SignedHeaders={signed_headers}, Signature={signature}")
    req_headers = {**all_headers, "Authorization": auth}
    return req_headers


def request(method, path, query_params=None, extra_headers=None):
    query = urllib.parse.urlencode(query_params or {})
    url = f"https://{HOST}{canonical_uri_for(path)}" + (f"?{query}" if query else "")
    req_headers = build_auth(method, path, query_params or {}, extra_headers or {})
    req = urllib.request.Request(url, method=method, headers={**req_headers, **(extra_headers or {})})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return resp.read()


def list_objects(prefix="", max_pages=None):
    keys, token = [], None
    pages = 0
    while True:
        pages += 1
        if max_pages and pages > max_pages:
            break
        params = {"list-type": "2", "max-keys": "1000"}
        if prefix:
            params["prefix"] = prefix
        if token:
            params["continuation-token"] = token
        xml = request("GET", f"/{BUCKET}", params)
        root = ET.fromstring(xml)
        ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        for c in root.findall(".//s3:Contents", ns):
            keys.append((c.findtext("s3:Key", namespaces=ns),
                         int(c.findtext("s3:Size", namespaces=ns) or 0)))
        token = root.findtext(".//s3:NextContinuationToken", namespaces=ns)
        if not token:
            break
    return keys


def get_object(key, dest, resume=False):
    import urllib.request as _ur
    offset = 0
    try:
        offset = os.path.getsize(dest)
    except OSError:
        offset = 0
    if resume and offset > 0:
        headers = {"Range": f"bytes={offset}-"}
        mode = "ab"
    else:
        headers = {}
        mode = "wb"
        offset = 0
    url = f"https://{HOST}{canonical_uri_for(f'/{BUCKET}/{key}')}"
    req_headers = build_auth("GET", f"/{BUCKET}/{key}", {}, headers)
    req = _ur.Request(url, headers={**req_headers, **headers})
    total = offset
    with _ur.urlopen(req, timeout=600) as resp:
        with open(dest, mode) as f:
            while True:
                chunk = resp.read(4 * 1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                total += len(chunk)
    return total


def list_dirs(prefix=""):
    dirs = []
    params = {"list-type": "2", "max-keys": "1000", "delimiter": "/"}
    if prefix:
        params["prefix"] = prefix
    while True:
        xml = request("GET", f"/{BUCKET}", params)
        root = ET.fromstring(xml)
        ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        for c in root.findall(".//s3:CommonPrefixes/s3:Prefix", ns):
            dirs.append(c.text)
        token = root.findtext(".//s3:NextContinuationToken", namespaces=ns)
        if not token:
            break
    return dirs


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "list"
    if cmd == "list":
        prefix = ""
        max_pages = None
        if "--prefix" in sys.argv:
            prefix = sys.argv[sys.argv.index("--prefix") + 1]
        if "--max-pages" in sys.argv:
            max_pages = int(sys.argv[sys.argv.index("--max-pages") + 1])
        keys = list_objects(prefix, max_pages)
        total = 0
        for k, sz in keys:
            total += sz
            print(f"{sz:>15,}  {k}")
        print(f"--- {len(keys)} objects, {total:,} bytes total")
    elif cmd == "dirs":
        prefix = ""
        if "--prefix" in sys.argv:
            prefix = sys.argv[sys.argv.index("--prefix") + 1]
        for d in list_dirs(prefix):
            print(d)
        print(f"--- dirs under '{prefix}'")
    elif cmd == "find":
        keyword = sys.argv[2]
        prefix = "final_somascan_smp/" if len(sys.argv) < 4 else sys.argv[3]
        max_pages = int(sys.argv[4]) if len(sys.argv) > 4 else None
        keys = list_objects(prefix, max_pages)
        hits = [(k, sz) for k, sz in keys if keyword.upper() in k.upper()]
        print(f"--- {len(hits)} hits for '{keyword}' in '{prefix}' (searched {len(keys)} keys)")
        for k, sz in hits:
            print(f"{sz:>15,}  {k}")
    elif cmd == "get":
        key, dest = sys.argv[2], sys.argv[3]
        resume = "--resume" in sys.argv
        n = get_object(key, dest, resume=resume)
        print(f"downloaded {key} -> {dest} ({n:,} bytes)")
    elif cmd == "sync":
        dest_dir = sys.argv[2]
        prefix = ""
        if "--prefix" in sys.argv:
            prefix = sys.argv[sys.argv.index("--prefix") + 1]
        os.makedirs(dest_dir, exist_ok=True)
        for key, sz in list_objects(prefix):
            dest = os.path.join(dest_dir, key.replace("/", os.sep))
            os.makedirs(os.path.dirname(dest) or dest_dir, exist_ok=True)
            if os.path.exists(dest) and os.path.getsize(dest) == sz:
                print(f"skip (exists) {key}")
                continue
            n = get_object(key, dest, resume=True)
            print(f"downloaded {key} ({n:,} bytes)")
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
