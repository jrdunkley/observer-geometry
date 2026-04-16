from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / "papers"
OUT = ROOT / "audit" / "tex_theorem_inventory.md"

ENV_RE = re.compile(
    r"\\begin\{(theorem|proposition|lemma|definition|corollary|remark|remarkplain|conjecture)\}"
    r"(?:\[(?P<title>[^\]]*)\])?"
    r"(?P<body>.*?)"
    r"\\end\{\1\}",
    re.DOTALL,
)
LABEL_RE = re.compile(r"\\label\{([^}]*)\}")


def clean_tex(text: str, limit: int = 1200) -> str:
    text = re.sub(r"%.*", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) > limit:
        text = text[: limit - 3].rstrip() + "..."
    return text


def main() -> None:
    lines: list[str] = ["# TeX Theorem Inventory", ""]
    for path in sorted(PAPERS.glob("*.tex")):
        text = path.read_text(encoding="utf-8", errors="replace")
        lines.append(f"## {path.name}")
        lines.append("")
        count = 0
        for match in ENV_RE.finditer(text):
            count += 1
            env = match.group(1)
            title = match.group("title") or ""
            body = match.group("body")
            label_match = LABEL_RE.search(body)
            label = label_match.group(1) if label_match else ""
            line_no = text.count("\n", 0, match.start()) + 1
            heading = f"- `{env}`"
            if title:
                heading += f" [{clean_tex(title, 160)}]"
            if label:
                heading += f" `{label}`"
            heading += f" (line {line_no})"
            lines.append(heading)
            lines.append(f"  - {clean_tex(body)}")
        if count == 0:
            lines.append("- no theorem-like environments found")
        lines.append("")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
