"""
release_app.py -- WSGI wrapper that adds GET /version and leaves the frozen core untouched.

Point your server entry (gunicorn / waitress / the WSGI file) at  release_app:app  instead of  server:app.
Every other path is passed straight to server.app, so the frozen engine and its SHA-256 are not changed.
"""
import json, os, time

import server  # the frozen core, imported unmodified

_HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_HERE, "RELEASE.json"), encoding="utf-8") as _f:
    RELEASE = json.load(_f)
_START = time.time()


def version_payload():
    short = getattr(server, "CORE_SHORT_SHA256", None)
    expected = getattr(server, "EXPECTED_BUILD_PREFIX", None)
    pinned = RELEASE["core_sha256"]
    return {
        "name": RELEASE["name"],
        "version": RELEASE["version"],
        "git_tag": RELEASE["git_tag"],
        "git_commit": os.environ.get("GIT_COMMIT") or RELEASE.get("git_commit"),
        "released": RELEASE.get("released"),
        "core_sha256_pinned": pinned,
        "core_short_sha256": short,
        "expected_build_prefix": expected,
        "build_verified": bool(short) and short == expected and pinned.startswith(short),
        "doi": RELEASE.get("doi"),
        "concept_doi": RELEASE.get("concept_doi"),
        "changelog": RELEASE.get("changelog"),
        "uptime_s": round(time.time() - _START, 1),
    }


def app(environ, start_response):
    if environ.get("PATH_INFO", "").rstrip("/") == "/version" and environ.get("REQUEST_METHOD") in ("GET", "HEAD"):
        body = json.dumps(version_payload(), indent=2).encode("utf-8")
        start_response("200 OK", [("Content-Type", "application/json"),
                                  ("Content-Length", str(len(body))),
                                  ("Cache-Control", "no-store")])
        return [b""] if environ.get("REQUEST_METHOD") == "HEAD" else [body]
    return server.app(environ, start_response)
