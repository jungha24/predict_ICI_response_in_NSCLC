from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Optional
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pandas as pd
import yaml

XML_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def log_step(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def ensure_outdir(path: str | Path) -> Path:
    outdir = Path(path)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def read_yaml(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_yaml(obj: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(obj, handle, sort_keys=False, allow_unicode=True)


def write_json(obj: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, ensure_ascii=False)


def normalize_id(value: object) -> str:
    text = str(value).strip().upper()
    text = re.sub(r"[_\s]+", "-", text)
    return re.sub(r"-+", "-", text)


def read_table(path: Optional[str], sheet_name: Optional[str | int] = 0) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    path_obj = Path(path)
    if not path_obj.exists():
        return None

    suffix = path_obj.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path_obj)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path_obj, sep="\t")
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path_obj)
    if suffix in {".xlsx", ".xlsm"}:
        try:
            return pd.read_excel(path_obj, sheet_name=sheet_name, engine="openpyxl")
        except ImportError:
            return read_xlsx_fallback(path_obj, sheet_name=sheet_name)
    raise ValueError(f"Unsupported file format: {path_obj}")


def read_xlsx_fallback(path: str | Path, sheet_name: str | int = 0) -> pd.DataFrame:
    workbook_path = Path(path)
    with ZipFile(workbook_path) as archive:
        shared_strings = _read_shared_strings(archive)
        sheet_target = _resolve_sheet_target(archive, sheet_name)
        rows = _read_sheet_rows(archive, sheet_target, shared_strings)

    if not rows:
        return pd.DataFrame()

    header = [str(cell) if cell is not None else "" for cell in rows[0]]
    header = make_unique_names([col if col else f"unnamed_{idx + 1}" for idx, col in enumerate(header)])
    body = rows[1:]
    if not body:
        return pd.DataFrame(columns=header)
    return pd.DataFrame(body, columns=header)


def make_unique_names(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for name in names:
        count = seen.get(name, 0)
        seen[name] = count + 1
        out.append(name if count == 0 else f"{name}_{count + 1}")
    return out


def sanitize_feature_name(name: str) -> str:
    text = re.sub(r"[^0-9A-Za-z_]+", "_", str(name).strip())
    return re.sub(r"_+", "_", text).strip("_")


def _read_shared_strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for si in root.findall(f"{XML_NS}si"):
        text_parts = [node.text or "" for node in si.iter(f"{XML_NS}t")]
        strings.append("".join(text_parts))
    return strings


def _resolve_sheet_target(archive: ZipFile, sheet_name: str | int) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}

    sheets = []
    for sheet in workbook.find(f"{XML_NS}sheets"):
        rel_id = sheet.attrib[f"{REL_NS}id"]
        sheets.append((sheet.attrib.get("name", ""), rel_map[rel_id]))

    if isinstance(sheet_name, str):
        matches = [target for name, target in sheets if name == sheet_name]
        if not matches:
            raise ValueError(f"Sheet '{sheet_name}' not found in {archive.filename}.")
        target = matches[0]
    else:
        sheet_index = int(sheet_name)
        if sheet_index < 0 or sheet_index >= len(sheets):
            raise ValueError(f"Sheet index {sheet_index} is out of range for {archive.filename}.")
        target = sheets[sheet_index][1]

    return target if target.startswith("xl/") else f"xl/{target}"


def _read_sheet_rows(archive: ZipFile, target: str, shared_strings: list[str]) -> list[list[object]]:
    root = ET.fromstring(archive.read(target))
    sheet_data = root.find(f"{XML_NS}sheetData")
    if sheet_data is None:
        return []

    rows: list[list[object]] = []
    max_width = 0
    for row in sheet_data.findall(f"{XML_NS}row"):
        row_values: dict[int, object] = {}
        for cell in row.findall(f"{XML_NS}c"):
            ref = cell.attrib.get("r", "")
            col_idx = _column_ref_to_index(ref)
            row_values[col_idx] = _parse_cell_value(cell, shared_strings)
            max_width = max(max_width, col_idx + 1)
        values = [None] * max_width
        for col_idx, value in row_values.items():
            if col_idx >= len(values):
                values.extend([None] * (col_idx + 1 - len(values)))
            values[col_idx] = value
        rows.append(values)

    for idx, row in enumerate(rows):
        if len(row) < max_width:
            rows[idx] = row + [None] * (max_width - len(row))
    return rows


def _column_ref_to_index(cell_ref: str) -> int:
    letters = "".join(char for char in cell_ref if char.isalpha()).upper()
    index = 0
    for char in letters:
        index = index * 26 + (ord(char) - ord("A") + 1)
    return max(index - 1, 0)


def _parse_cell_value(cell: ET.Element, shared_strings: list[str]) -> object:
    cell_type = cell.attrib.get("t")
    value_node = cell.find(f"{XML_NS}v")
    inline_node = cell.find(f"{XML_NS}is")

    if inline_node is not None:
        text_parts = [node.text or "" for node in inline_node.iter(f"{XML_NS}t")]
        return "".join(text_parts)

    if value_node is None:
        return None

    raw_value = value_node.text
    if raw_value is None:
        return None
    if cell_type == "s":
        return shared_strings[int(raw_value)]
    if cell_type == "b":
        return int(raw_value)

    try:
        numeric_value = float(raw_value)
        if numeric_value.is_integer():
            return int(numeric_value)
        return numeric_value
    except ValueError:
        return raw_value
