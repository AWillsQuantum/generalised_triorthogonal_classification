"""Compact finite domains of pointed stabiliser supports."""

from collections import Counter
from contextlib import contextmanager
import struct

from affine_graph import rank
from compact_to_native import BASE_RECORD, MANIFEST_HEADER
from space_codec import read_exact
from verify_large_quotient_domain import pointed_support


DOMAIN_HEADER = struct.Struct("<8sIQQI")
POINTING = struct.Struct("<IBi")
MAGIC = b"UTPDOM1\0"


def write_header(stream, source_count, count, length_limit):
    if (type(source_count) is not int or not 0 < source_count < 1 << 32
            or type(count) is not int or not 0 <= count < 1 << 64
            or type(length_limit) is not int or not 0 < length_limit <= 64):
        raise ValueError("Invalid finite pointing domain")
    stream.write(DOMAIN_HEADER.pack(MAGIC, 1, source_count, count, length_limit))


def read_header(stream):
    magic, version, sources, count, limit = DOMAIN_HEADER.unpack(read_exact(stream, DOMAIN_HEADER.size))
    if magic != MAGIC or version != 1 or not 0 < sources < 1 << 32 or not 0 < limit <= 64:
        raise ValueError("Invalid pointing-domain header")
    return dict(source_count=sources, count=count, maximum_protocol_length=limit)


def read_pointings(stream, header):
    for _ in range(header["count"]):
        source, parity, origin = POINTING.unpack(read_exact(stream, POINTING.size))
        if source >= header["source_count"] or parity not in range(5) or (parity < 3 and origin < 0):
            raise ValueError("Invalid pointing-domain record")
        yield source, parity, origin
    if stream.read(1):
        raise ValueError("Trailing pointing-domain data")


@contextmanager
def open_domain(path):
    import zstandard
    with path.open("rb") as raw:
        with zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
            yield decoded


def source_key(source):
    row = source["space"]
    return row["c"], row["m"], row["index"]


def source_points(source, supports):
    key = source_key(source)
    points = supports[key]
    frame = source.get("coordinate_frame")
    if frame is None:
        return points
    m = key[1]
    if (len(frame) != m+1 or any(type(x) is not int or not 0 <= x < 1 << m for x in frame)
            or rank(frame[1:]) != m):
        raise ValueError("Invalid source coordinate frame")
    def image(point):
        value = frame[0]
        for i in range(m):
            if point >> i & 1:
                value ^= frame[i+1]
        return value
    return tuple(sorted(map(image, points)))


def reconstruct(sources, supports, record):
    source, parity, origin = record
    row = sources[source]
    key = source_key(row)
    return pointed_support(source_points(row, supports), key[1], parity, origin)


def write_native(stream, header, sources, supports, records):
    """Reconstruct the native input, retaining its finite record order."""
    count = header["count"]
    if len(sources) != header["source_count"]:
        raise ValueError("Source domain and pointing header disagree")
    stream.write(MANIFEST_HEADER.pack(b"UTSPTS1\0", 1, len(sources), count,
                                     header["maximum_protocol_length"], 1))
    for index, source in enumerate(sources):
        if source["index"] != index:
            raise ValueError("Nonconsecutive source identities")
        c, m, i = source_key(source)
        for value in (f"c{c}_m{m}_i{i}", f"c{c}_m{m}"):
            raw = value.encode("ascii")
            stream.write(struct.pack("<H", len(raw)))
            stream.write(raw)
    seen = 0
    for seen, record in enumerate(records, 1):
        if seen > count:
            raise ValueError("Too many pointing records")
        source, parity, origin = record
        c, m, _ = source_key(sources[source])
        points, h = reconstruct(sources, supports, record)
        if len(points) > header["maximum_protocol_length"]:
            raise ValueError("Pointing exceeds the protocol-length scope")
        stream.write(BASE_RECORD.pack(source, c, m, h, parity, len(points), origin))
        stream.write(struct.pack(f"<{len(points)}I", *points))
    if seen != count:
        raise ValueError("Too few pointing records")


def verify_records(stream, sources, supports):
    header = read_header(stream)
    if header["source_count"] != len(sources):
        raise ValueError("Inconsistent source count")
    counts, parities, lengths = Counter(), Counter(), Counter()
    for record in read_pointings(stream, header):
        points, _ = reconstruct(sources, supports, record)
        if not 0 < len(points) <= header["maximum_protocol_length"]:
            raise ValueError("Pointing outside protocol-length scope")
        counts[record[0]] += 1
        parities[record[1]] += 1
        lengths[len(points)] += 1
    if set(counts) != set(range(len(sources))):
        raise ValueError("Source without a pointed representative")
    return dict(**header, pointing_counts_by_source=dict(counts),
                parity_counts=dict(parities), protocol_length_counts=dict(lengths))


def write_base_native(stream, sources, supports, length_limit):
    count = len(sources)
    stream.write(MANIFEST_HEADER.pack(b"UTSPTS1\0", 1, count, count, length_limit, 3))
    for index, source in enumerate(sources):
        if source["index"] != index:
            raise ValueError("Nonconsecutive source indices")
        c, m, i = source_key(source)
        for value in (f"c{c}_m{m}_i{i}", f"c{c}_m{m}"):
            raw = value.encode("ascii")
            stream.write(struct.pack("<H", len(raw)))
            stream.write(raw)
    for index, source in enumerate(sources):
        c, m, _ = source_key(source)
        points = source_points(source, supports)
        stream.write(BASE_RECORD.pack(index, c, m, m, 0, c, 0))
        stream.write(struct.pack(f"<{c}I", *points))


def read_native_header(stream):
    magic, version, count, records, limit, flags = MANIFEST_HEADER.unpack(read_exact(stream, MANIFEST_HEADER.size))
    if magic != b"UTSPTS1\0" or version not in (1, 2) or flags not in range(4):
        raise ValueError("Invalid native support header")
    names = []
    for _ in range(count):
        pair = []
        for _ in range(2):
            size, = struct.unpack("<H", read_exact(stream, 2))
            pair.append(read_exact(stream, size).decode("ascii"))
        names.append(tuple(pair))
    return dict(version=version, source_count=count, count=records,
                maximum_protocol_length=limit, flags=flags, names=names)


def read_native_records(stream, header):
    for index in range(header["count"]):
        source, c, m, h, parity, n, origin = BASE_RECORD.unpack(read_exact(stream, BASE_RECORD.size))
        original_index = index
        if header["version"] == 2:
            original_index, = struct.unpack("<Q", read_exact(stream, 8))
        points = struct.unpack(f"<{n}I", read_exact(stream, 4*n))
        if source >= header["source_count"] or parity not in range(5):
            raise ValueError("Invalid native source or pointing case")
        yield dict(source_index=source, c=c, m=m, ambient_dimension=h,
                   parity_case=parity, origin=origin, points=points,
                   original_record_index=original_index)
    if stream.read(1):
        raise ValueError("Trailing native support data")
