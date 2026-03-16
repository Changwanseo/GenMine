import socket
import time
import urllib.request
import urllib.error
import ssl
import logging


NCBI_HOST = "eutils.ncbi.nlm.nih.gov"
NCBI_TEST_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
NCBI_SEARCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=nucleotide&term=fungi&retmax=1"
)


def check_dns():
    logging.info("[1/5] DNS Resolution...")
    try:
        ip = socket.gethostbyname(NCBI_HOST)
        logging.info(f"  OK: {NCBI_HOST} -> {ip}")
        return True
    except socket.gaierror as e:
        logging.error(f"  FAIL: Cannot resolve {NCBI_HOST}: {e}")
        logging.error("  -> DNS may be blocked or misconfigured")
        return False


def check_tcp_connection():
    logging.info("[2/5] TCP Connection (port 443)...")
    try:
        start = time.time()
        sock = socket.create_connection((NCBI_HOST, 443), timeout=10)
        elapsed = time.time() - start
        sock.close()
        logging.info(f"  OK: Connected in {elapsed:.2f}s")
        if elapsed > 2:
            logging.warning(f"  WARNING: High latency ({elapsed:.2f}s) - may cause timeouts")
        return True
    except socket.timeout:
        logging.error("  FAIL: Connection timed out (10s)")
        logging.error("  -> Port 443 may be blocked by firewall")
        return False
    except OSError as e:
        logging.error(f"  FAIL: {e}")
        return False


def check_ssl_handshake():
    logging.info("[3/5] SSL/TLS Handshake...")
    try:
        context = ssl.create_default_context()
        start = time.time()
        with socket.create_connection((NCBI_HOST, 443), timeout=10) as sock:
            with context.wrap_socket(sock, server_hostname=NCBI_HOST) as ssock:
                elapsed = time.time() - start
                proto = ssock.version()
                logging.info(f"  OK: {proto} handshake in {elapsed:.2f}s")
                return True
    except ssl.SSLError as e:
        logging.error(f"  FAIL: SSL error: {e}")
        logging.error("  -> University firewall may be intercepting HTTPS (DPI/MITM)")
        return False
    except Exception as e:
        logging.error(f"  FAIL: {e}")
        return False


def check_http_request(api_key=None):
    logging.info("[4/5] HTTPS Request (einfo endpoint, up to 3 attempts)...")
    for attempt in range(3):
        try:
            start = time.time()
            url = NCBI_TEST_URL
            if api_key:
                url += f"?api_key={api_key}"
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "GenMine/NetworkCheck")
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                elapsed = time.time() - start
                status = resp.status
                if attempt > 0:
                    logging.info(f"  OK on attempt {attempt+1}: HTTP {status}, {len(data)} bytes in {elapsed:.2f}s")
                else:
                    logging.info(f"  OK: HTTP {status}, {len(data)} bytes in {elapsed:.2f}s")
                if elapsed > 3:
                    logging.warning(f"  WARNING: Response took {elapsed:.2f}s (expected <1s) - slow connection to NCBI")
                return True
        except urllib.error.HTTPError as e:
            logging.warning(f"  Attempt {attempt+1}/3: HTTP {e.code}: {e.reason}")
            if e.code == 500 and attempt < 2:
                time.sleep(2)
                continue
            if e.code == 403:
                logging.error("  -> Access forbidden - possible firewall/proxy block")
            elif e.code == 429:
                logging.error("  -> Rate limited by NCBI")
            elif e.code == 500:
                logging.error("  -> NCBI server error on all attempts - server may be overloaded")
            return False
        except urllib.error.URLError as e:
            logging.error(f"  FAIL: {e.reason}")
            logging.error("  -> Network unreachable or blocked by proxy/firewall")
            return False
        except Exception as e:
            logging.error(f"  FAIL: {e}")
            return False
    return False


def check_repeated_requests(api_key=None):
    logging.info("[5/5] Sustained Connection Test (5 rapid requests)...")
    results = []
    for i in range(5):
        try:
            start = time.time()
            url = NCBI_SEARCH_URL
            if api_key:
                url += f"&api_key={api_key}"
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "GenMine/NetworkCheck")
            with urllib.request.urlopen(req, timeout=15) as resp:
                resp.read()
                elapsed = time.time() - start
                results.append(("OK", elapsed))
                logging.info(f"  Request {i+1}/5: OK ({elapsed:.2f}s)")
        except Exception as e:
            elapsed = time.time() - start
            results.append(("FAIL", elapsed))
            logging.error(f"  Request {i+1}/5: FAIL after {elapsed:.2f}s - {e}")
        time.sleep(0.4)

    ok_count = sum(1 for r, _ in results if r == "OK")
    if ok_count == 5:
        times = [t for r, t in results if r == "OK"]
        avg = sum(times) / len(times)
        logging.info(f"  OK: All 5 requests succeeded (avg {avg:.2f}s)")
    elif ok_count > 0:
        logging.warning(f"  PARTIAL: {ok_count}/5 succeeded - connection is UNSTABLE")
        logging.warning("  -> Intermittent drops suggest firewall throttling or packet loss")
    else:
        logging.error("  FAIL: All 5 requests failed")
        logging.error("  -> NCBI is completely unreachable from this network")
    return ok_count


def run_network_check(api_key=None):
    logging.info("=" * 60)
    logging.info("GenMine Network Diagnostics")
    logging.info(f"Target: {NCBI_HOST}")
    if api_key:
        logging.info(f"API Key: {api_key[:4]}...{api_key[-4:]}")
    else:
        logging.info("API Key: not provided (rate limit: 3 req/s)")
    logging.info("=" * 60)

    dns_ok = check_dns()
    if not dns_ok:
        logging.error("")
        logging.error("DIAGNOSIS: DNS resolution failed.")
        logging.error("  - Check your internet connection")
        logging.error("  - Try: nslookup eutils.ncbi.nlm.nih.gov")
        logging.error("  - University DNS may be blocking external lookups")
        logging.error("  - Workaround: use 8.8.8.8 or 1.1.1.1 as DNS server")
        return False

    tcp_ok = check_tcp_connection()
    if not tcp_ok:
        logging.error("")
        logging.error("DIAGNOSIS: TCP connection to port 443 blocked.")
        logging.error("  - University firewall is likely blocking outbound HTTPS")
        logging.error("  - Contact IT to whitelist eutils.ncbi.nlm.nih.gov")
        logging.error("  - Workaround: use VPN or mobile hotspot")
        return False

    ssl_ok = check_ssl_handshake()
    if not ssl_ok:
        logging.error("")
        logging.error("DIAGNOSIS: SSL/TLS handshake failed.")
        logging.error("  - University may be running HTTPS inspection (DPI/MITM proxy)")
        logging.error("  - This intercepts encrypted traffic and can break API connections")
        logging.error("  - Contact IT or use VPN to bypass inspection")
        return False

    http_ok = check_http_request(api_key=api_key)
    if not http_ok:
        logging.error("")
        logging.error("DIAGNOSIS: HTTP request failed despite TCP+SSL working.")
        logging.error("  - Possible HTTP-level proxy filtering")
        logging.error("  - Check if you need to authenticate with a captive portal")
        logging.error("  - Try opening https://eutils.ncbi.nlm.nih.gov in a browser")
        return False

    ok_count = check_repeated_requests(api_key=api_key)

    logging.info("")
    logging.info("=" * 60)
    if ok_count == 5:
        logging.info("RESULT: Network connection to NCBI is HEALTHY")
        logging.info("GenMine should work normally on this network.")
    elif ok_count >= 3:
        logging.warning("RESULT: Network connection is UNSTABLE")
        logging.warning("GenMine may fail intermittently. Recommendations:")
        logging.warning("  - Use VPN or mobile hotspot for reliable downloads")
        logging.warning("  - Run GenMine during off-peak hours (late night/early morning)")
        logging.warning("  - Contact university IT about connection stability to NCBI")
    elif ok_count > 0:
        logging.warning("RESULT: Network connection is VERY UNSTABLE")
        logging.warning("GenMine will likely fail. Recommendations:")
        logging.warning("  - Use VPN or mobile hotspot")
        logging.warning("  - Contact university IT to whitelist eutils.ncbi.nlm.nih.gov")
    else:
        logging.error("RESULT: NCBI is UNREACHABLE")
        logging.error("GenMine cannot function on this network without a workaround.")
    logging.info("=" * 60)
    return ok_count == 5
