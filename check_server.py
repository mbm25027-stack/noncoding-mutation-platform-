#!/usr/bin/env python3
"""
Server health checker (v3.3) - detects DOWN vs FROZEN servers.

States:
  OK       - responded in time with the expected status (and version, if set)
  SLOW     - responded, but slower than --slow seconds
  FROZEN   - TCP port accepts connections, but HTTP reply times out (hung process)
  DOWN     - cannot connect at all (process dead / port closed / network issue)
  BAD      - responded, but wrong status code or wrong version

Usage examples:
  python check_server.py http://localhost:8000/health
  python check_server.py http://localhost:8000/health --version 3.3 --interval 30
  python check_server.py http://localhost:8000/health --interval 30 \
      --restart-cmd "systemctl restart myserver" --fail-limit 3

Only uses the Python standard library (no installs needed).
"""
import argparse
import json
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from urllib.parse import urlparse


def tcp_open(host: str, port: int, timeout: float) -> bool:
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except OSError:
        return False


def check_once(url: str, timeout: float, slow: float, expect_status: int,
               expect_version: str | None) -> tuple[str, str]:
    parsed = urlparse(url)
    host = parsed.hostname or "localhost"
    port = parsed.port or (443 if parsed.scheme == "https" else 80)

    start = time.monotonic()
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "health-check/3.3"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read(65536).decode("utf-8", errors="replace")
            status = resp.status
    except urllib.error.HTTPError as e:
        return "BAD", f"HTTP {e.code} (expected {expect_status})"
    except (socket.timeout, TimeoutError):
        # Timed out waiting. Is the port still accepting connections?
        if tcp_open(host, port, min(timeout, 3)):
            return "FROZEN", f"port {port} open but no HTTP reply within {timeout}s"
        return "DOWN", f"timed out and port {port} not accepting connections"
    except urllib.error.URLError as e:
        reason = e.reason
        if isinstance(reason, (socket.timeout, TimeoutError)):
            if tcp_open(host, port, min(timeout, 3)):
                return "FROZEN", f"port {port} open but no HTTP reply within {timeout}s"
            return "DOWN", "timed out and port not accepting connections"
        return "DOWN", f"connection failed: {reason}"
    except Exception as e:  # noqa: BLE001
        return "DOWN", f"error: {e}"

    elapsed = time.monotonic() - start

    if status != expect_status:
        return "BAD", f"HTTP {status} (expected {expect_status})"

    if expect_version:
        found = None
        try:
            found = str(json.loads(body).get("version"))
        except (ValueError, AttributeError):
            pass
        if found is None:
            found = body.strip()[:50]
        if expect_version not in found:
            return "BAD", f"version mismatch: got '{found}', expected '{expect_version}'"

    if elapsed > slow:
        return "SLOW", f"OK but took {elapsed:.2f}s (limit {slow}s)"
    return "OK", f"responded in {elapsed:.2f}s"


def main() -> int:
    p = argparse.ArgumentParser(description="Server health checker v3.3")
    p.add_argument("url", help="health URL, e.g. http://localhost:8000/health")
    p.add_argument("--timeout", type=float, default=5.0, help="seconds before a request counts as hung")
    p.add_argument("--slow", type=float, default=2.0, help="seconds above which reply is SLOW")
    p.add_argument("--status", type=int, default=200, help="expected HTTP status")
    p.add_argument("--version", default=None, help="expected version string, e.g. 3.3")
    p.add_argument("--interval", type=float, default=0, help="repeat every N seconds (0 = run once)")
    p.add_argument("--fail-limit", type=int, default=3, help="consecutive failures before restart")
    p.add_argument("--restart-cmd", default=None, help="shell command to run when server is FROZEN/DOWN")
    args = p.parse_args()

    failures = 0
    while True:
        state, detail = check_once(args.url, args.timeout, args.slow, args.status, args.version)
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {state:<6} {detail}", flush=True)

        if state in ("FROZEN", "DOWN"):
            failures += 1
        else:
            failures = 0

        if args.restart_cmd and failures >= args.fail_limit:
            print(f"  -> {failures} failures in a row, running restart command", flush=True)
            subprocess.run(args.restart_cmd, shell=True, check=False)
            failures = 0
            time.sleep(10)  # give the server time to boot

        if args.interval <= 0:
            return 0 if state in ("OK", "SLOW") else 1
        time.sleep(args.interval)


if __name__ == "__main__":
    sys.exit(main())