from basalt.qc_backend import CheckM2Backend, CheckMBackend


def test_checkm2_command_quotes_paths():
    command = CheckM2Backend().build_cmd(4, "input folder", "output folder")
    assert "-i 'input folder'" in command
    assert "-o 'output folder'" in command


def test_checkm_command_quotes_paths():
    command = CheckMBackend().build_cmd(4, "input folder", "output folder")
    assert "'input folder' 'output folder'" in command
