from scripts.utils import sanitize_seqid


def test_sanitize_seqid_basic() -> None:
    assert sanitize_seqid("NM_001178324.2") == "NM_001178324.2"
    assert sanitize_seqid("my seq/id") == "my_seq_id"
