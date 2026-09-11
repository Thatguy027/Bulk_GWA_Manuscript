#!/usr/bin/env python3
"""Ask a 9.5 GB remote PRIDE result file which SID-2 peptides were observed.

  python3 scripts/query_sid2_psms_pride.py

WHY THIS EXISTS
METHODS claims that the peptide carrying the N94-C95-T96 sequon was never
identified in Tan et al. 2024's intestinal proteome, and that SID-2's whole
confident detection there is two cytoplasmic-tail peptides. That claim needs
the PSM table, which is not in the paper's supplement -- only protein-level
counts are. It is in the deposited Proteome Discoverer result file, which is
9.53 GB.

Downloading 9.53 GB to read about twenty rows is not necessary. A .msf file is
a bare SQLite database and PRIDE answers HTTP range requests
(Accept-Ranges: bytes), so this script implements just enough of the SQLite
file format -- page reads, varint decoding, record decoding, table and index
b-tree descent, overflow-page assembly -- to query the file where it sits. The
whole run costs roughly 570 range requests and 2.3 MB, 0.02% of the file.

THE QUERY PATH, which is why it is cheap
  sqlite_master (page 1)          -> schema, table and index root pages
  TargetProteins                  -> scan for accession G5EEV9
  IX_TargetProteinsTargetPsms_... -> index descent on the protein id, giving
                                     that protein's (WorkflowID, PeptideID)
                                     links; 4 page reads for 21 links
  IX_TargetPsms_Unique_...        -> index descent per link, giving a rowid
  TargetPsms                      -> rowid descent for each PSM row

Only the TargetProteins scan is linear, and it stops at the first match.

WHAT IT FINDS, and the assertions below pin it
Five of the 21 SID-2 PSMs clear a 1% Percolator q value. All five are from
intestine runs and none from gonad runs, which reproduces the paper's tissue
call at the level of single spectra. They are two peptides, ITVPDAMR (239-246)
and TNGLYGYDNNNSSR (225-238), both cytoplasmic: 22 residues, 7.1% coverage,
matching the 2 peptides / 5 PSMs / 7% coverage in the deposited report. The
ectodomain is unsampled -- the only ectodomain peptides present at all are
SKPFQLGFASATLNR (q = 0.21) and NVTFILEVTTDTK (q = 0.61), both from gonad runs.
NCTFTANYTGYFTPDPK, which spans N94-C95-T96, is absent in any form.

Needs network. Reads the deposited SID-2 sequence for residue numbering, so it
runs from a clone; nothing else in the figure pipeline depends on it.
---------------------------------------------------------------------------"""
import json, os, re, struct, sys

try:
    import requests
except ImportError:
    sys.exit("needs `requests`")

URL = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/05/PXD047792/"
       "James_20220518-G-I_chimerys20230616.msf")
FASTA = "supplemental_data/structure/sid2_ortholog_sequences.fa"
ACC = "G5EEV9"
ECD = (21, 193)
SEQUON_PEPTIDE = "NCTFTANYTGYFTPDPK"      # C. elegans 94-110
Q_CUT = 0.01


# ---------------------------------------------------------------- SQLite ----
class Pager:
    """Read-only SQLite pager backed by HTTP range requests."""

    def __init__(self, url):
        self.url, self.s = url, requests.Session()
        self.cache, self.bytes_read, self.reqs = {}, 0, 0
        hdr = self._range(0, 99)
        if hdr[:15] != b"SQLite format 3":
            raise SystemExit("remote file is not a SQLite database")
        ps = struct.unpack(">H", hdr[16:18])[0]
        self.page_size = 65536 if ps == 1 else ps
        self.n_pages = struct.unpack(">I", hdr[28:32])[0]
        self.usable = self.page_size - hdr[20]

    def _range(self, a, b):
        r = self.s.get(self.url, headers={"Range": f"bytes={a}-{b}"}, timeout=120)
        r.raise_for_status()
        self.bytes_read += len(r.content)
        self.reqs += 1
        return r.content

    def page(self, n):
        if n not in self.cache:
            off = (n - 1) * self.page_size
            self.cache[n] = self._range(off, off + self.page_size - 1)
        return self.cache[n]


def varint(buf, i):
    v = 0
    for k in range(9):
        b = buf[i + k]
        if k == 8:
            return (v << 8) | b, i + 9
        v = (v << 7) | (b & 0x7F)
        if not b & 0x80:
            return v, i + k + 1


_FIX = {0: (0, "null"), 1: (1, "int"), 2: (2, "int"), 3: (3, "int"), 4: (4, "int"),
        5: (6, "int"), 6: (8, "int"), 7: (8, "float"), 8: (0, "zero"), 9: (0, "one")}


def rec_decode(pay):
    hlen, i = varint(pay, 0)
    types, j = [], i
    while j < hlen:
        t, j = varint(pay, j)
        types.append(t)
    out, k = [], hlen
    for t in types:
        if t in _FIX:
            n, kind = _FIX[t]
            raw, k = pay[k:k + n], k + n
            out.append({"null": None, "zero": 0, "one": 1}.get(kind) if kind in
                       ("null", "zero", "one") else
                       (int.from_bytes(raw, "big", signed=True) if kind == "int"
                        else struct.unpack(">d", raw)[0]))
        elif t >= 12 and t % 2 == 0:
            n = (t - 12) // 2
            out.append(pay[k:k + n]); k += n
        elif t >= 13:
            n = (t - 13) // 2
            out.append(pay[k:k + n].decode("utf-8", "replace")); k += n
        else:
            out.append(None)
    return out


def _assemble(pg, p, off, plen, xmax):
    """Payload, following overflow pages. xmax is the page type's local max."""
    U = pg.usable
    if plen <= xmax:
        return p[off:off + plen]
    M = ((U - 12) * 32 // 255) - 23
    K = M + ((plen - M) % (U - 4))
    local = K if K <= xmax else M
    buf = bytearray(p[off:off + local])
    nxt = struct.unpack(">I", p[off + local:off + local + 4])[0]
    while nxt and len(buf) < plen:
        op = pg.page(nxt)
        nxt = struct.unpack(">I", op[0:4])[0]
        buf += op[4:4 + min(U - 4, plen - len(buf))]
    return bytes(buf)


def _cellptrs(p, base, hdr):
    n = struct.unpack(">H", p[base + 3:base + 5])[0]
    return [struct.unpack(">H", p[base + hdr + 2 * i:base + hdr + 2 * i + 2])[0]
            for i in range(n)]


def table_scan(pg, root):
    """Yield (rowid, values) over a table b-tree, depth first."""
    p = pg.page(root)
    base = 100 if root == 1 else 0
    typ = p[base]
    if typ == 0x05:
        for q in _cellptrs(p, base, 12):
            yield from table_scan(pg, struct.unpack(">I", p[q:q + 4])[0])
        right = struct.unpack(">I", p[base + 8:base + 12])[0]
        if right:
            yield from table_scan(pg, right)
    elif typ == 0x0D:
        for q in _cellptrs(p, base, 8):
            plen, i = varint(p, q)
            rowid, i = varint(p, i)
            yield rowid, rec_decode(_assemble(pg, p, i, plen, pg.usable - 35))


def index_find(pg, root, key):
    """Yield index records whose leading columns equal `key`, pruning subtrees."""
    key = tuple(key) if isinstance(key, (tuple, list)) else (key,)
    xmax = ((pg.usable - 12) * 64 // 255) - 23

    def cells(p, base, typ):
        hdr = 12 if typ == 0x02 else 8
        out = []
        for q in _cellptrs(p, base, hdr):
            if typ == 0x02:
                child = struct.unpack(">I", p[q:q + 4])[0]
                plen, i = varint(p, q + 4)
            else:
                child, (plen, i) = None, varint(p, q)
            out.append((child, rec_decode(_assemble(pg, p, i, plen, xmax))))
        return out

    def rec(pageno):
        p = pg.page(pageno)
        base = 100 if pageno == 1 else 0
        typ = p[base]
        cs = cells(p, base, typ)
        if typ == 0x0A:
            for _, r in cs:
                if tuple(r[:len(key)]) == key:
                    yield r
            return
        prev = None
        for child, r in cs:
            ck = tuple(r[:len(key)])
            if (prev is None or prev <= key) and key <= ck:
                yield from rec(child)
            if ck == key:
                yield r
            prev = ck
        if prev is None or key >= prev:
            right = struct.unpack(">I", p[base + 8:base + 12])[0]
            if right:
                yield from rec(right)

    yield from rec(root)


def row_by_rowid(pg, root, rowid):
    pageno = root
    while True:
        p = pg.page(pageno)
        base = 100 if pageno == 1 else 0
        typ = p[base]
        if typ == 0x05:
            nxt = None
            for q in _cellptrs(p, base, 12):
                k, _ = varint(p, q + 4)
                if rowid <= k:
                    nxt = struct.unpack(">I", p[q:q + 4])[0]
                    break
            pageno = nxt or struct.unpack(">I", p[base + 8:base + 12])[0]
        else:
            for q in _cellptrs(p, base, 8):
                plen, i = varint(p, q)
                rid, i = varint(p, i)
                if rid == rowid:
                    return rec_decode(_assemble(pg, p, i, plen, pg.usable - 35))
            return None


# ------------------------------------------------------------------ query ----
def sid2_sequence():
    seq, take = [], False
    for line in open(FASTA):
        if line.startswith(">"):
            take = line.startswith(">" + ACC)
            if not take and seq:
                break
            continue
        if take:
            seq.append(line.strip())
    return "".join(seq)


def main():
    if not os.path.exists(FASTA):
        sys.exit(f"run from the repository root; {FASTA} not found")
    S = sid2_sequence()
    print(f"SID-2 {ACC}: {len(S)} aa, ectodomain {ECD[0]}-{ECD[1]}")

    pg = Pager(URL)
    print(f"remote: {pg.n_pages * pg.page_size / 1e9:.2f} GB, "
          f"page size {pg.page_size}")

    schema, cols = {}, {}
    for _, v in table_scan(pg, 1):
        if v[0] == "table":
            schema[v[1]] = v[3]
            if v[1] == "TargetPsms":
                cols = {c.strip().split()[0].strip('[]"'): i for i, c in
                        enumerate(v[4][v[4].find("(") + 1:].split(","))}
        elif v[0] == "index":
            schema[v[1]] = v[3]
    ix_prot = next(k for k in schema if k.startswith(
        "IX_TargetProteinsTargetPsms_Unique_"))
    ix_psm = "IX_TargetPsms_Unique_WorkflowID_PeptideID"

    prot_id = None
    for _, v in table_scan(pg, schema["TargetProteins"]):
        if v[3] and ACC in str(v[3]):
            prot_id = v[0]
            break
    if prot_id is None:
        sys.exit(f"{ACC} not in TargetProteins")
    print(f"{ACC} UniqueSequenceID = {prot_id}")

    links = [(r[1], r[2]) for r in index_find(pg, schema[ix_prot], prot_id)]
    psms = []
    for wf, pid in links:
        hit = next(iter(index_find(pg, schema[ix_psm], (wf, pid))), None)
        if hit:
            r = row_by_rowid(pg, schema["TargetPsms"], hit[-1])
            if r:
                psms.append(r)
    print(f"{len(psms)} SID-2 PSMs, read in {pg.reqs} range requests "
          f"({pg.bytes_read / 1e6:.2f} MB, "
          f"{100 * pg.bytes_read / (pg.n_pages * pg.page_size):.3f}% of the file)\n")

    def tissue(r):
        tag = str(r[cols["SpectrumFileName"]]).split("Aur_1h_")[-1]
        return "intestine" if tag.startswith("Ix") else "gonad"

    rows = []
    for r in psms:
        seq = r[cols["Sequence"]]
        start = S.find(seq) + 1
        rows.append(dict(seq=seq, start=start, end=start + len(seq) - 1,
                         q=r[cols["PercolatorqValue"]], tissue=tissue(r),
                         mods=r[cols["Modifications"]],
                         region=("cytoplasmic" if start > ECD[1] else
                                 "ectodomain" if start + len(seq) - 1 <= ECD[1]
                                 else "spans TM")))
    print(f"{'peptide':20s} {'residues':>10s} {'q':>10s}  {'tissue':9s} region")
    for d in sorted(rows, key=lambda d: (d["seq"], d["q"])):
        flag = " *" if d["q"] < Q_CUT else ""
        print(f"{d['seq']:20s} {d['start']:4d}-{d['end']:<5d} {d['q']:10.4g}  "
              f"{d['tissue']:9s} {d['region']}{flag}")

    conf = [d for d in rows if d["q"] < Q_CUT]
    peps = {d["seq"] for d in conf}
    covered = {i for d in conf for i in range(d["start"], d["end"] + 1)}
    print(f"\n* clears q < {Q_CUT}: {len(conf)} PSMs, {len(peps)} peptides, "
          f"{len(covered)} residues = {100 * len(covered) / len(S):.1f}% coverage")
    print(f"  tissues of those PSMs: "
          f"{sorted({d['tissue'] for d in conf})}")
    print(f"  regions: {sorted({d['region'] for d in conf})}")
    print(f"  {SEQUON_PEPTIDE} (N94-C95-T96) observed: "
          f"{any(d['seq'] == SEQUON_PEPTIDE for d in rows)}")

    # the claims METHODS makes, so a change in the deposit fails loudly
    assert len(conf) == 5, "the confident PSM count changed"
    assert peps == {"ITVPDAMR", "TNGLYGYDNNNSSR"}, "the confident peptides changed"
    assert {d["tissue"] for d in conf} == {"intestine"}, \
        "a confident SID-2 PSM now comes from a gonad run"
    assert {d["region"] for d in conf} == {"cytoplasmic"}, \
        "a confident SID-2 peptide is no longer cytoplasmic"
    assert round(100 * len(covered) / len(S), 1) == 7.1, "coverage changed"
    assert not any(d["seq"] == SEQUON_PEPTIDE for d in rows), \
        "the N94 sequon peptide now appears -- METHODS must be rewritten"
    print("\nall assertions hold")


if __name__ == "__main__":
    main()
