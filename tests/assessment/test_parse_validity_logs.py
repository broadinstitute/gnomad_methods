"""Tests for the parse_validity_logs module."""

from gnomad.assessment.parse_validity_logs import generate_html_report, parse_log_file


def _write_log(tmp_path, text):
    """Write log text to a file and return its path."""
    path = str(tmp_path / "validity.log")
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return path


def _read_file(path):
    """Read and return the full contents of a file."""
    with open(path, encoding="utf-8") as f:
        return f.read()


class TestParseLogFile:
    """Test the parse_log_file function."""

    def test_categorizes_levels_and_maps_function(self, tmp_path):
        """Test log-level/keyword categorization and function-name mapping."""
        log = "\n".join(
            [
                "INFO (gnomad.x.check_sex_chr_metrics 1): passed the check",
                "INFO (gnomad.x.check_missingness 2): this check failed",
                "INFO (gnomad.x.some_other_func 3): just an update",
                "WARNING (gnomad.x.check_missingness 4): heads up",
                "ERROR (gnomad.x.check_missingness 5): broke",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        # (validity_check, category, source, message, table)
        categories = [(p[0], p[1]) for p in parsed]
        assert categories == [
            ("XY check", "pass"),
            ("missingness", "fail"),
            ("some_other_func", "info"),  # unmapped name passes through
            ("missingness", "warn"),
            ("missingness", "fail"),
        ]
        assert parsed[0][2] == "gnomad.x.check_sex_chr_metrics 1"

    def test_associates_table_with_preceding_log(self, tmp_path):
        """Test that an ASCII table is captured with its preceding log entry."""
        log = "\n".join(
            [
                "INFO (gnomad.x.check_missingness 10): found issues",
                "+-------+",
                "| locus |",
                "+-------+",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        assert len(parsed) == 1
        message, table = parsed[0][3], parsed[0][4]
        assert message == "found issues"
        assert "| locus |" in table

    def test_accumulates_multiline_message(self, tmp_path):
        """Test that continuation lines are joined into one entry's message."""
        log = "\n".join(
            [
                "INFO (gnomad.x.main 1): first line",
                "second line",
                "third line",
                "INFO (gnomad.x.main 2): next entry",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        assert len(parsed) == 2
        assert parsed[0][3] == "first line\nsecond line\nthird line"
        assert parsed[1][3] == "next entry"

    def test_fail_keyword_uses_word_boundaries(self, tmp_path):
        """Test that 'fail' matches only as a word and 'passed' takes precedence."""
        log = "\n".join(
            [
                "INFO (gnomad.x.main 1): no failures detected",
                "INFO (gnomad.x.main 2): check failed but later passed",
                "INFO (gnomad.x.main 3): the check failed",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        categories = [p[1] for p in parsed]
        # "failures" is not the word "fail"/"failed"; "passed" wins over "failed".
        assert categories == ["info", "pass", "fail"]

    def test_empty_file_returns_empty_list(self, tmp_path):
        """Test that an empty log file yields no parsed entries."""
        assert parse_log_file(_write_log(tmp_path, "")) == []

    def test_ignores_lines_before_first_entry(self, tmp_path):
        """Test that non-matching lines before any log entry are dropped."""
        log = "\n".join(
            [
                "preamble noise that is not a log line",
                "another stray line",
                "INFO (gnomad.x.main 1): real entry",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        assert len(parsed) == 1
        assert parsed[0][3] == "real entry"

    def test_table_associated_with_non_final_entry(self, tmp_path):
        """Test that a table is captured for an entry followed by another entry.

        This exercises the in-loop flush path (a new match closing out a prior
        entry that had a table), distinct from the end-of-file flush.
        """
        log = "\n".join(
            [
                "INFO (gnomad.x.check_missingness 1): found issues",
                "+-------+",
                "| locus |",
                "+-------+",
                "INFO (gnomad.x.main 2): all done",
            ]
        )
        parsed = parse_log_file(_write_log(tmp_path, log))

        assert len(parsed) == 2
        # First entry keeps its message and captures the table.
        assert parsed[0][3] == "found issues"
        assert "| locus |" in parsed[0][4]
        # The trailing entry has no table.
        assert parsed[1][3] == "all done"
        assert parsed[1][4] == ""


class TestGenerateHtmlReport:
    """Test the generate_html_report function."""

    def test_escapes_content_and_builds_filters(self, tmp_path):
        """Test HTML escaping of messages and per-column filter options."""
        parsed_logs = [
            (
                "missingness",
                "fail",
                "gnomad.x.check_missingness 1",
                "broke <b> & died",
                "",
            ),
            ("XY check", "pass", "gnomad.x.check_sex_chr_metrics 2", "all good", ""),
        ]
        out = str(tmp_path / "report.html")
        generate_html_report(parsed_logs, out)
        html_out = _read_file(out)

        # Message content is HTML-escaped, not injected raw.
        assert "broke &lt;b&gt; &amp; died" in html_out
        assert "broke <b>" not in html_out
        # Filter dropdowns list each validity check and status.
        assert '<option value="missingness">missingness</option>' in html_out
        assert '<option value="XY check">XY check</option>' in html_out
        assert '<option value="fail">FAIL</option>' in html_out
        assert '<option value="pass">PASS</option>' in html_out
        # Category drives the row CSS class.
        assert '<tr class="fail">' in html_out

    def test_view_table_button_only_when_table_present(self, tmp_path):
        """Test that a "View Table" toggle is rendered only for entries with tables."""
        parsed_logs = [
            ("missingness", "fail", "src 1", "has table", "+---+\n| x |\n+---+"),
            ("XY check", "pass", "src 2", "no table", ""),
        ]
        out = str(tmp_path / "report.html")
        generate_html_report(parsed_logs, out)
        html_out = _read_file(out)

        # Exactly one toggle button/table div, for the single entry with a table.
        assert html_out.count('class="toggle-btn"') == 1
        assert '<div id="table_0" class="hidden-table">' in html_out
        assert 'id="table_1"' not in html_out
