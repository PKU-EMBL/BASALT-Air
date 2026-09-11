from pathlib import Path

from basalt import cli
from basalt.core.cleanup import cleanup


def test_short_read_parser_accepts_relative_directory_paths():
    assert cli._parse_paired_reads("reads/r1.fq,reads/r2.fq") == {
        "1": ["reads/r1.fq", "reads/r2.fq"]
    }


def test_short_read_parser_keeps_legacy_group_separator():
    assert cli._parse_paired_reads("r1.fq,r2.fq/s1.fq,s2.fq") == {
        "1": ["r1.fq", "r2.fq"],
        "2": ["s1.fq", "s2.fq"],
    }


def test_input_paths_are_anchored_before_workdir_change(tmp_path, monkeypatch):
    source = tmp_path / "data"
    source.mkdir()
    for name in ("assembly.fa", "r1.fq", "r2.fq"):
        (source / name).write_text("input")

    monkeypatch.setattr(cli, "assembly_list", ["data/assembly.fa"], raising=False)
    monkeypatch.setattr(cli, "lr_list", [], raising=False)
    monkeypatch.setattr(cli, "hifi_list", [], raising=False)
    monkeypatch.setattr(cli, "hic_list", [], raising=False)
    monkeypatch.setattr(cli, "coverage_list", [], raising=False)
    monkeypatch.setattr(cli, "datasets", {"1": ["data/r1.fq", "data/r2.fq"]}, raising=False)
    monkeypatch.setattr(cli, "data_feeding_folder", [], raising=False)
    monkeypatch.setattr(cli, "binsets_list", [], raising=False)
    monkeypatch.setattr(cli, "refinement_binset", "", raising=False)

    cli._resolve_input_paths(str(tmp_path))
    workdir = tmp_path / "work"
    workdir.mkdir()
    monkeypatch.chdir(workdir)
    cli._normalize_paths_to_cwd()

    assert Path("assembly.fa").resolve() == (source / "assembly.fa").resolve()
    assert Path("r1.fq").resolve() == (source / "r1.fq").resolve()
    assert Path("r2.fq").resolve() == (source / "r2.fq").resolve()


def test_cleanup_preserves_pipeline_log(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path("Basalt_log.txt").write_text("completed")
    cleanup([])
    assert Path("Basalt_log.txt").read_text() == "completed"
