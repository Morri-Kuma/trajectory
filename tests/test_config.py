from src.utils import load_config, resolve_path

def test_load_config_test_profile():
    cfg = load_config("config.yaml")
    assert cfg["mode"] == "test"
    assert cfg["active"] is cfg["test"]
    assert len(cfg["_config_hash"]) >= 6

def test_server_profile_present_and_scvi_on():
    cfg = load_config("config.yaml")
    assert cfg["server"]["scvi"]["use_scvi"] is True
    assert cfg["annotation"]["reference_label_key"] == "author_cluster_label"

def test_resolve_path():
    cfg = {"project": {"root": "."}}
    assert resolve_path(cfg, "a/b").endswith("a" + __import__("os").sep + "b")
