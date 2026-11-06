"""Deterministic local text retrieval for an AIMD Analysis Kit workspace."""

from __future__ import annotations

import json
import math
import re
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1
DEFAULT_INDEX_NAME = "index.json"
DEFAULT_CHUNK_LINES = 80
DEFAULT_OVERLAP_LINES = 20
DEFAULT_MAX_FILE_BYTES = 2_000_000
DEFAULT_EXTENSIONS = frozenset(
    {
        ".csv",
        ".json",
        ".md",
        ".ps1",
        ".py",
        ".rst",
        ".sh",
        ".toml",
        ".txt",
        ".yml",
        ".yaml",
    }
)
DEFAULT_FILENAMES = frozenset({".gitignore", "CHANGELOG", "LICENSE", "Makefile", "README"})
DEFAULT_EXCLUDED_DIRS = frozenset(
    {
        ".git",
        ".github",
        ".pytest_cache",
        ".retrieval",
        ".ruff_cache",
        ".venv",
        "__pycache__",
        "_smoke",
        "_staging",
        "build",
        "dist",
        "node_modules",
        "results",
        "venv",
    }
)

# Keep identifiers and paths intact while making Chinese text searchable one character at a time.
TOKEN_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_./:+-]*|[\u4e00-\u9fff]")


class RetrievalError(ValueError):
    """Raised for malformed indexes or invalid retrieval options."""


@dataclass(frozen=True)
class RetrievalResult:
    path: str
    start_line: int
    end_line: int
    score: float
    snippet: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "start_line": self.start_line,
            "end_line": self.end_line,
            "score": round(self.score, 6),
            "snippet": self.snippet,
        }


def tokenize(text: str) -> list[str]:
    """Return case-folded tokens while keeping technical identifiers searchable."""

    return [match.group(0).casefold() for match in TOKEN_RE.finditer(text)]


def _normalise_extensions(extensions: set[str] | frozenset[str] | list[str] | None) -> set[str]:
    values = DEFAULT_EXTENSIONS if extensions is None else extensions
    normalised = {value.casefold() if value.startswith(".") else f".{value.casefold()}" for value in values}
    if not normalised:
        raise RetrievalError("at least one file extension must be selected")
    return normalised


def _is_candidate(path: Path, root: Path, extensions: set[str]) -> bool:
    relative_parts = path.relative_to(root).parts
    if any(part in DEFAULT_EXCLUDED_DIRS for part in relative_parts[:-1]):
        return False
    return path.name in DEFAULT_FILENAMES or path.suffix.casefold() in extensions


def _read_text(path: Path, max_file_bytes: int) -> str | None:
    if path.stat().st_size > max_file_bytes:
        return None
    raw = path.read_bytes()
    if b"\x00" in raw:
        return None
    return raw.decode("utf-8", errors="replace")


def _chunk_file(text: str, chunk_lines: int, overlap_lines: int) -> list[tuple[int, int, str]]:
    lines = text.splitlines()
    if not lines:
        return []
    step = chunk_lines - overlap_lines
    chunks: list[tuple[int, int, str]] = []
    for start in range(0, len(lines), step):
        stop = min(start + chunk_lines, len(lines))
        content = "\n".join(lines[start:stop]).strip()
        if content:
            chunks.append((start + 1, stop, content))
        if stop == len(lines):
            break
    return chunks


def build_index(
    root: Path,
    index_path: Path,
    *,
    extensions: set[str] | frozenset[str] | list[str] | None = None,
    chunk_lines: int = DEFAULT_CHUNK_LINES,
    overlap_lines: int = DEFAULT_OVERLAP_LINES,
    max_file_bytes: int = DEFAULT_MAX_FILE_BYTES,
) -> dict[str, Any]:
    """Build and write a portable workspace index, returning its metadata."""

    root = root.expanduser().resolve()
    index_path = index_path.expanduser().resolve()
    if not root.is_dir():
        raise RetrievalError(f"retrieval root is not a directory: {root}")
    if chunk_lines < 1 or overlap_lines < 0 or overlap_lines >= chunk_lines:
        raise RetrievalError("chunk-lines must be positive and overlap-lines must be smaller")
    if max_file_bytes < 1:
        raise RetrievalError("max-file-bytes must be positive")
    selected_extensions = _normalise_extensions(extensions)

    chunks: list[dict[str, Any]] = []
    document_frequency: Counter[str] = Counter()
    document_count = 0
    skipped_count = 0
    for path in sorted(root.rglob("*"), key=lambda item: item.as_posix().casefold()):
        if not path.is_file() or path.is_symlink() or not _is_candidate(path, root, selected_extensions):
            continue
        try:
            text = _read_text(path, max_file_bytes)
        except OSError:
            skipped_count += 1
            continue
        if text is None:
            skipped_count += 1
            continue
        file_chunks = _chunk_file(text, chunk_lines, overlap_lines)
        if not file_chunks:
            continue
        document_count += 1
        relative_path = path.relative_to(root).as_posix()
        for start_line, end_line, content in file_chunks:
            frequencies = Counter(tokenize(content))
            if not frequencies:
                continue
            document_frequency.update(frequencies.keys())
            chunks.append(
                {
                    "path": relative_path,
                    "start_line": start_line,
                    "end_line": end_line,
                    "text": content,
                    "term_frequency": dict(frequencies),
                    "length": sum(frequencies.values()),
                }
            )

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root_name": root.name or root.drive,
        "document_count": document_count,
        "chunk_count": len(chunks),
        "skipped_count": skipped_count,
        "chunk_lines": chunk_lines,
        "overlap_lines": overlap_lines,
        "max_file_bytes": max_file_bytes,
        "extensions": sorted(selected_extensions),
        "document_frequency": dict(document_frequency),
        "chunks": chunks,
    }
    index_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=index_path.parent, delete=False) as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temporary_path = Path(handle.name)
    temporary_path.replace(index_path)
    return metadata


def load_index(index_path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(index_path.expanduser().read_text(encoding="utf-8"))
    except OSError as exc:
        raise RetrievalError(f"cannot read retrieval index: {index_path}") from exc
    except json.JSONDecodeError as exc:
        raise RetrievalError(f"invalid retrieval index JSON: {index_path}") from exc
    if payload.get("schema_version") != SCHEMA_VERSION or not isinstance(payload.get("chunks"), list):
        raise RetrievalError(f"unsupported retrieval index: {index_path}")
    return payload


def _snippet(text: str, query_tokens: set[str], max_lines: int = 8) -> str:
    lines = text.splitlines()
    matching = [index for index, line in enumerate(lines) if query_tokens.intersection(tokenize(line))]
    if not matching:
        return "\n".join(lines[:max_lines])
    center = matching[0]
    start = max(0, center - max_lines // 3)
    return "\n".join(lines[start : start + max_lines])


def search_index(index: dict[str, Any], query: str, *, top_k: int = 10) -> list[RetrievalResult]:
    """Search an index using a compact BM25-style score with path and phrase boosts."""

    if top_k < 1:
        raise RetrievalError("top-k must be positive")
    query_tokens = tokenize(query)
    if not query_tokens:
        raise RetrievalError("query must contain searchable text")
    chunks = index["chunks"]
    if not chunks:
        return []
    query_terms = set(query_tokens)
    frequencies = index.get("document_frequency", {})
    total_chunks = len(chunks)
    average_length = sum(int(chunk.get("length", 0)) for chunk in chunks) / total_chunks
    if average_length <= 0:
        return []
    results: list[RetrievalResult] = []
    for chunk in chunks:
        term_frequency = chunk.get("term_frequency", {})
        length = max(1, int(chunk.get("length", 0)))
        score = 0.0
        for term in query_terms:
            term_count = int(term_frequency.get(term, 0))
            if not term_count:
                continue
            doc_frequency = int(frequencies.get(term, 0))
            idf = math.log(1.0 + (total_chunks - doc_frequency + 0.5) / (doc_frequency + 0.5))
            denominator = term_count + 1.2 * (0.25 + 0.75 * length / average_length)
            score += idf * term_count * 2.2 / denominator
        if score <= 0:
            continue
        path_tokens = set(tokenize(chunk["path"]))
        score += 0.25 * len(query_terms.intersection(path_tokens))
        if query.casefold() in chunk["text"].casefold():
            score += 0.5
        results.append(
            RetrievalResult(
                path=chunk["path"],
                start_line=int(chunk["start_line"]),
                end_line=int(chunk["end_line"]),
                score=score,
                snippet=_snippet(chunk["text"], query_terms),
            )
        )
    results.sort(key=lambda item: (-item.score, item.path, item.start_line))
    return results[:top_k]


def format_results(results: list[RetrievalResult], output_format: str = "text") -> str:
    if output_format == "json":
        return json.dumps([result.as_dict() for result in results], ensure_ascii=False, indent=2)
    if not results:
        return "No matching workspace text found."
    rendered: list[str] = []
    for rank, result in enumerate(results, start=1):
        rendered.append(f"{rank}. {result.path}:{result.start_line}-{result.end_line} score={result.score:.3f}")
        rendered.append(result.snippet)
        rendered.append("")
    return "\n".join(rendered).rstrip()
