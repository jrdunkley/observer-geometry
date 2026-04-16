from __future__ import annotations

import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from .engine import fetch_coinbase_hourly, run_observer_gate, synthetic_hourly_candles


ROOT = Path(__file__).resolve().parent
STATIC = ROOT / "static"


class Handler(BaseHTTPRequestHandler):
    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, status: int, payload: dict) -> None:
        self._send(status, json.dumps(payload, indent=2).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path in {"/", "/index.html"}:
            body = (STATIC / "index.html").read_bytes()
            self._send(200, body, "text/html; charset=utf-8")
            return
        if parsed.path == "/api/run":
            try:
                query = parse_qs(parsed.query)
                limit = int(query.get("limit", ["300"])[0])
                train_window = int(query.get("train_window", ["96"])[0])
                abstain_z = float(query.get("abstain_z", ["0.15"])[0])
                mode = query.get("mode", ["live"])[0]
                candles = synthetic_hourly_candles(limit) if mode == "synthetic" else fetch_coinbase_hourly(limit=limit)
                result = run_observer_gate(
                    candles,
                    train_window=train_window,
                    abstain_z=abstain_z,
                    source="synthetic" if mode == "synthetic" else "coinbase",
                )
                result["mode"] = "synthetic_offline" if mode == "synthetic" else "live_coinbase"
                self._send_json(200, result)
            except Exception as exc:  # pragma: no cover - user-facing server guard
                self._send_json(500, {"error": type(exc).__name__, "message": str(exc)})
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def log_message(self, format: str, *args) -> None:  # noqa: A002
        print("%s - %s" % (self.address_string(), format % args))


def main() -> None:
    server = ThreadingHTTPServer(("127.0.0.1", 8765), Handler)
    print("BTC observer gate UI: http://127.0.0.1:8765")
    print("Use Ctrl+C to stop.")
    server.serve_forever()


if __name__ == "__main__":
    main()
