#!/usr/bin/env python3
"""Convert a standard MATPOWER case file to GridSolve-C's text format."""

import argparse
import re
from pathlib import Path


def matrix(text, name):
    match = re.search(r"mpc\." + name + r"\s*=\s*\[(.*?)\]\s*;", text, re.S)
    if not match:
        raise ValueError(f"missing mpc.{name} matrix")
    rows = []
    for line in match.group(1).splitlines():
        line = line.split("%", 1)[0].strip().rstrip(";")
        if line:
            rows.append([float(value) for value in line.split()])
    return rows


def scalar(text, name, default=None):
    match = re.search(r"mpc\." + name + r"\s*=\s*([-+0-9.eE]+)\s*;", text)
    if not match:
        if default is not None:
            return default
        raise ValueError(f"missing mpc.{name}")
    return float(match.group(1))


def convert(source, destination):
    text = Path(source).read_text(encoding="utf-8")
    base = scalar(text, "baseMVA", 100.0)
    buses = matrix(text, "bus")
    generators = matrix(text, "gen")
    branches = matrix(text, "branch")
    gens = {}
    for row in generators:
        if row[7] == 0:
            continue
        bus = int(row[0])
        entry = gens.setdefault(bus, [0.0, 0.0, 0.0, 0.0])
        entry[0] += row[1]
        entry[1] += row[2]
        entry[2] += row[3]
        entry[3] += row[4]
    lines = [f"BASEMVA {base:.12g}"]
    for row in buses:
        bus_id, bus_type = int(row[0]), int(row[1])
        pg, qg, qmax, qmin = gens.get(bus_id, [0.0, 0.0, 1e9, -1e9])
        pd, qd, gs, bs = row[2:6]
        vm, va, vmax, vmin = row[7], row[8], row[11], row[12]
        lines.append("BUS " + " ".join(f"{value:.12g}" for value in
            [bus_id, bus_type, vm, va, pd, qd, pg, qg, gs/base, bs/base, vmin, vmax, qmin, qmax]))
    for row in branches:
        ratio = row[8] if row[8] != 0 else 1.0
        lines.append("BRANCH " + " ".join(f"{value:.12g}" for value in
            [int(row[0]), int(row[1]), row[2], row[3], row[4], ratio, row[9], row[5], int(row[10])]))
    Path(destination).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="MATPOWER case .m file")
    parser.add_argument("output", help="GridSolve-C .grid output path")
    args = parser.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()