import pytest

from basalt.shell import CommandError, run_cmd


def test_run_cmd_supports_argument_lists():
    result = run_cmd(["python", "-c", "print('ok')"], capture=True)
    assert result.stdout.strip() == "ok"


def test_run_cmd_reports_missing_executable():
    with pytest.raises(CommandError, match="unable to execute"):
        run_cmd(["basalt-command-that-does-not-exist"])
