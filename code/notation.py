"""Current resource notation with lossless support for certified legacy data.

The historical enumeration kernels and frozen evidence use S, q and h.
Current public witnesses use N, k and r, with N = k + r.
"""

from copy import deepcopy

LEGACY_SCHEMA = "triorthogonal-protocol-witnesses-v1"
CURRENT_SCHEMA = "triorthogonal-protocol-witnesses-v2"
DEFINITION = "N is the total number of matrix rows; k is the number of logical rows; r=N-k is the number of stabiliser rows."


def current_catalogue(data):
    """Change names and gate typography, never matrices or output identities."""
    if data.get("schema") == CURRENT_SCHEMA:
        legacy_catalogue(data)
        return deepcopy(data)
    if data.get("schema") != LEGACY_SCHEMA:
        raise ValueError("Unsupported witness schema")
    result = deepcopy(data)
    result["schema"] = CURRENT_SCHEMA
    result["pareto_objectives"] = ["n", "N"]
    result["footprint_definition"] = DEFINITION
    result["notation_update"] = dict(
        schema="resource-notation-update-v1", mapping={"S": "N", "q": "k", "h": "r"},
        previous_footprint_definition=data["footprint_definition"],
        stable_identifiers="Output IDs, matrix entries, source paths and hashes are unchanged.",
        frozen_evidence="Original files retain their original notation; see NOTATION.md.")
    for row in result["outputs"]:
        row["k"] = row.pop("q")
        for factor in row.get("independent_factors", []):
            factor["k"] = factor.pop("q")
        if "factorised_representative_latex" in row:
            row["factorised_representative_latex"] = row["factorised_representative_latex"].replace(
                r"\mathsf{", r"\mathrm{")
    for row in result["protocols"]:
        row["k"] = row.pop("q")
        row["N"] = row.pop("S")
        row["r"] = row["N"] - row["k"]
    return result


def legacy_catalogue(data):
    """Validate the notation layer and return the historical mathematical schema."""
    if data.get("schema") == LEGACY_SCHEMA:
        return deepcopy(data)
    if data.get("schema") != CURRENT_SCHEMA:
        raise ValueError("Unsupported witness schema")
    if data.get("pareto_objectives") != ["n", "N"] or data.get("footprint_definition") != DEFINITION:
        raise ValueError("Incorrect resource definition")
    result = deepcopy(data)
    metadata = result.pop("notation_update")
    if metadata.get("mapping") != {"S": "N", "q": "k", "h": "r"}:
        raise ValueError("Incorrect notation mapping")
    result["schema"] = LEGACY_SCHEMA
    result["pareto_objectives"] = ["n", "S"]
    result["footprint_definition"] = metadata["previous_footprint_definition"]
    for row in result["outputs"]:
        if "q" in row:
            raise ValueError("Conflicting old and new logical dimensions")
        row["q"] = row.pop("k")
        for factor in row.get("independent_factors", []):
            if "q" in factor:
                raise ValueError("Conflicting factor dimensions")
            factor["q"] = factor.pop("k")
        if "factorised_representative_latex" in row:
            row["factorised_representative_latex"] = row["factorised_representative_latex"].replace(
                r"\mathrm{", r"\mathsf{")
    for row in result["protocols"]:
        if "q" in row or "S" in row:
            raise ValueError("Conflicting old and new resources")
        k, N, r = row["k"], row["N"], row.pop("r")
        if any(type(v) is not int for v in (k, N, r)) or not 0 < k <= N or r != N-k:
            raise ValueError("Inconsistent N, k, r")
        row["q"] = row.pop("k")
        row["S"] = row.pop("N")
    if current_catalogue(result) != data:
        raise ValueError("Unexpected changes in the notation layer")
    return result


def resource_parameters(rows, k):
    """Return resource labels for an independent-row protocol matrix."""
    if type(k) is not int or not 0 < k <= len(rows) or not rows or not rows[0]:
        raise ValueError("Invalid logical dimension or empty matrix")
    if any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("Ragged matrix")
    return dict(n=len(rows[0]), N=len(rows), k=k, r=len(rows)-k)
