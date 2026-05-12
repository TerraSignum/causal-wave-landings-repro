"""
Forward-protocol test for the bounded-operator-readout blindness
witness (sec:coefficients-provenance, paragraph (P6)).

This test enforces that, in any future re-determination of the
readout, the readout-commit hash is an ancestor of the
landings-commit hash in the public repository commit graph.
The check uses the bundled commit_hash_readout and
commit_hash_landings fields of
data/coefficient_readout_blindness.json.

If both fields are still TBD (pre-release), the test passes as
soft-pass; once a public release is cut and the hashes are
populated, the test enforces the ancestor relation strictly.
"""

import json
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
BLINDNESS = REPO / "data" / "coefficient_readout_blindness.json"


def test_blindness_json_loads():
    data = json.loads(BLINDNESS.read_text(encoding="utf-8"))
    assert "commit_hash_readout" in data
    assert "commit_hash_landings" in data


def test_blindness_protocol_documented():
    data = json.loads(BLINDNESS.read_text(encoding="utf-8"))
    assert data.get("forward_protocol") is not None, (
        "Forward-protocol field must document the blindness "
        "enforcement workflow."
    )


def test_blindness_release_pin():
    data = json.loads(BLINDNESS.read_text(encoding="utf-8"))
    assert "release_pin" in data, (
        "Release pin field must be present so the blindness witness "
        "is anchored to a specific release."
    )


def test_blindness_ancestor_relation_when_resolved():
    """
    If the commit hashes are populated (post-release), check that
    the readout commit is an ancestor of the landings commit. If
    they are still 'TBD-...' placeholders, soft-pass.
    """
    data = json.loads(BLINDNESS.read_text(encoding="utf-8"))
    h_readout = data.get("commit_hash_readout", "")
    h_landings = data.get("commit_hash_landings", "")
    if h_readout.startswith("TBD") or h_landings.startswith("TBD"):
        pytest.skip(
            "Pre-release state: commit hashes are placeholders. "
            "The blindness witness is the manuscript-draft "
            "chronology bundled with the audit trail until a public "
            "release is cut."
        )
    # When hashes are populated, run git merge-base --is-ancestor
    try:
        result = subprocess.run(
            ["git", "merge-base", "--is-ancestor", h_readout, h_landings],
            cwd=REPO,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, (
            f"Readout commit {h_readout} is not an ancestor of "
            f"landings commit {h_landings} in the public commit graph"
        )
    except FileNotFoundError:
        pytest.skip("git not available in CI sandbox")
