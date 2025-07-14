#!/usr/bin/env python3
"""Strip unwanted metadata keys from CWL files and ensure `sd:version: 100` on workflows.

Usage:  strip_cwl_metadata.py file1.cwl file2.cwl ...
"""
import sys, pathlib, ruamel.yaml as ry
BAD_ROOT_KEYS = {"namespaces", "schemas", "s:mainEntity", "s:name", "s:alternateName", "s:downloadUrl", "s:codeRepository", "s:license", "s:isPartOf", "s:creator", "s:about"}
BAD_PORT_KEYS = {"format"}

yaml = ry.YAML()
yaml.preserve_quotes = True
yaml.indent(mapping=2, sequence=4, offset=2)

def clean_cwl(path: pathlib.Path):
    data = yaml.load(path.read_text())
    # 1. Remove unwanted root-level keys
    for key in list(data.keys()):
        if key in BAD_ROOT_KEYS:
            del data[key]
    # 2. Remove unwanted keys from inputs and outputs
    for section in ("inputs", "outputs"):
        if section in data and isinstance(data[section], dict):
            for port in data[section].values():
                if isinstance(port, dict):
                    for bad in BAD_PORT_KEYS:
                        if bad in port:
                            del port[bad]
    # 3. Ensure sd:version for workflows
    if path.parts[0] == "workflows" and "sd:version" not in data:
        data["sd:version"] = 100
    # Write back
    with path.open('w', encoding='utf-8') as fh:
        yaml.dump(data, fh)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: strip_cwl_metadata.py <cwl files>", file=sys.stderr)
        sys.exit(1)
    for f in sys.argv[1:]:
        clean_cwl(pathlib.Path(f))
