"""Script to parse logs from validity check output."""

import re

import hail as hl


def parse_log_file(log_file):
    """Parse a log file and categorizes messages for formatting, extracting function names and sources.

    :param log_file: Path to the log file containing python logging output (such as INFO and WARNING statements).
    :return: List of tuples containing logger information describing the validity check, status category, source, message, and associated table if relevant.
    """
    parsed_logs = []
    log_pattern = re.compile(
        r"^(INFO|WARNING|ERROR) \(([^)]+)\.(\w+) (\d+)\): (.*)"
    )  # Extract log level, module, function, line number, and message.

    function_mapping = {
        "validate_config": "general info",
        "validate_ht_fields": "general info",
        "main": "general info",
        "validate_federated_data": "general info",
        "summarize_variants": "variant summary",
        "sum_group_callstats": "group summations",
        "check_missingness": "missingness",
        "compute_missingness": "missingness",
        "make_group_sum_expr_dict": "group summations",
        "check_sex_chr_metrics": "XY check",
        "check_raw_and_adj_callstats": "raw/adj check",
    }

    with hl.hadoop_open(log_file, "r") as f:
        current_message = ""
        current_table = []
        current_metadata = None  # Stores log metadata before a table appears.

        for line in f:
            match = log_pattern.match(line)

            if match:
                # Store previous log if it had a table.
                if current_message and current_metadata:
                    # Remove ASCII table lines from the message part but keep the full
                    #  table separately.
                    cleaned_message = re.sub(r"\+\-.*\+$", "", current_message).strip()
                    parsed_logs.append(
                        (
                            *current_metadata,
                            cleaned_message,
                            "\n".join(current_table).strip(),
                        )
                    )

                # Start new log message.
                log_level, module, function_name, line_number, message = match.groups()
                source = f"{module}.{function_name} {line_number}"
                validity_check = function_mapping.get(function_name, function_name)

                # Determine the category.
                message_lower = message.lower()
                log_categories = {
                    "INFO": lambda msg: (
                        "pass"
                        if "passed" in msg
                        else (
                            "fail"
                            if any(word in msg for word in ["failed", "fail"])
                            else "info"
                        )
                    ),
                    "WARNING": "warn",
                    "ERROR": "fail",
                }
                category = log_categories.get(log_level, "info")
                category = category(message_lower) if callable(category) else category

                # Reset tracking.
                current_message = message
                current_table = []
                current_metadata = (validity_check, category, source)

            elif "+----" in line or "| locus" in line:  # Table start detection.
                current_table.append(line.strip())  # Add table row.
            elif current_table:  # If already collecting a table.
                current_table.append(line.strip())

        # Store last log if it had a table.
        if current_message and current_metadata:
            cleaned_message = re.sub(r"\+\-.*\+$", "", current_message).strip()
            parsed_logs.append(
                (*current_metadata, cleaned_message, "\n".join(current_table).strip())
            )

    return parsed_logs


def generate_html_report(parsed_logs, output_file):
    """Generate an HTML report with sortable and filterable columns, with expandable tables for results."""
    html_template = """
    <html>
    <head>
        <style>
            body { font-family: Arial, sans-serif; }
            table { width: 100%; border-collapse: collapse; }
            th, td { border: 1px solid black; padding: 8px; text-align: left; vertical-align: top; }
            th { background-color: #f2f2f2; cursor: pointer; }
            .pass { color: 0d1cb6; }
            .fail { color: #D42736; font-weight: bold;}
            .warn { color: #DAA520; }
            .info { color: black; }
            .hidden-table { display: none; }
            .toggle-btn {
                cursor: pointer;
                color: blue;
                text-decoration: underline;
                float: right;  /* Moves "View Table" to the right */
                margin-left: 15px; /* Adds spacing between message and button */
                font-weight: normal;
            }
            .checkbox-container {
                display: inline-block;
                margin-left: 20px;
            }
            pre {
                text-align: left;
                white-space: pre-wrap;
                font-family: monospace;
                background: #f8f8f8;
                padding: 10px;
                border-radius: 5px;
                overflow-x: auto;
                width: 100%;
                display: block;
                color: black; /* Ensure table text is plain black */
                font-weight: normal; /* Ensure no bolding */
            }
        </style>
        <script>
            function toggleTable(id) {
                var tableDiv = document.getElementById(id);
                tableDiv.style.display = (tableDiv.style.display === "none" || tableDiv.style.display === "") ? "block" : "none";
            }

            function toggleAllTables() {
                var tables = document.getElementsByClassName("hidden-table");
                var checkbox = document.getElementById("toggleAll");
                var showAll = checkbox.checked;

                for (var i = 0; i < tables.length; i++) {
                    tables[i].style.display = showAll ? "block" : "none";
                }
            }

            function sortTable(n) {
                var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
                table = document.getElementById("logTable");
                switching = true;
                dir = "asc";
                while (switching) {
                    switching = false;
                    rows = table.rows;
                    for (i = 1; i < (rows.length - 1); i++) {
                        shouldSwitch = false;
                        x = rows[i].getElementsByTagName("TD")[n].innerHTML.toLowerCase();
                        y = rows[i + 1].getElementsByTagName("TD")[n].innerHTML.toLowerCase();
                        if ((dir == "asc" && x > y) || (dir == "desc" && x < y)) {
                            shouldSwitch = true;
                            break;
                        }
                    }
                    if (shouldSwitch) {
                        rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                        switching = true;
                        switchcount++;
                    } else {
                        if (switchcount === 0 && dir === "asc") {
                            dir = "desc";
                            switching = true;
                        }
                    }
                }
            }

            function filterTable() {
                var validityFilter = document.getElementById("functionFilter").value.toLowerCase();
                var statusFilter = document.getElementById("statusFilter").value.toLowerCase();
                var table, tr, i;
                table = document.getElementById("logTable");
                tr = table.getElementsByTagName("tr");

                for (i = 1; i < tr.length; i++) {
                    var validityCheck = tr[i].getElementsByTagName("td")[0].innerHTML.toLowerCase();
                    var status = tr[i].getElementsByTagName("td")[1].innerHTML.toLowerCase();

                    if ((validityFilter === "all" || validityCheck === validityFilter) &&
                        (statusFilter === "all" || status === statusFilter)) {
                        tr[i].style.display = "";
                    } else {
                        tr[i].style.display = "none";
                    }
                }
            }
        </script>
    </head>
    <body>
        <h2>Log Report</h2>
        <label for="functionFilter">Filter by Validity Check:</label>
        <select id="functionFilter" onchange="filterTable()">
            <option value="all">All</option>
    """

    validity_checks = set()
    statuses = set()

    for validity_check, category, source, message, table in parsed_logs:
        validity_checks.add(validity_check)
        statuses.add(category)

    for validity_check in sorted(validity_checks):
        html_template += f'<option value="{validity_check}">{validity_check}</option>'

    html_template += """
        </select>
        <label for="statusFilter">Filter by Status:</label>
        <select id="statusFilter" onchange="filterTable()">
            <option value="all">All</option>
    """

    for status in sorted(statuses):
        html_template += f'<option value="{status}">{status.upper()}</option>'

    html_template += """
        </select>

        <!-- Checkbox to show/hide all tables -->
        <span class="checkbox-container">
            <input type="checkbox" id="toggleAll" onclick="toggleAllTables()">
            <label for="toggleAll">Show All Tables</label>
        </span>

        <table id="logTable">
            <tr>
                <th onclick="sortTable(0)">Validity Check</th>
                <th onclick="sortTable(1)">Status</th>
                <th onclick="sortTable(2)">Source</th>
                <th onclick="sortTable(3)">Message</th>
            </tr>
    """

    for i, (validity_check, category, source, message, table) in enumerate(parsed_logs):
        table_id = f"table_{i}"
        table_button = (
            f'<span class="toggle-btn" onclick="toggleTable(\'{table_id}\')">View Table</span>'
            if table
            else ""
        )

        html_template += (
            f'<tr class="{category}">'
            f"<td>{validity_check}</td>"
            f'<td class="{category}">{category.upper()}</td>'
            f"<td>{source}</td>"
            f"<td>{message} {table_button}"
        )

        if table:
            html_template += (
                f'<div id="{table_id}" class="hidden-table"><pre>{table}</pre></div>'
            )

        html_template += "</td></tr>"

    html_template += "</table></body></html>"

    with hl.hadoop_open(output_file, "w") as f:
        f.write(html_template)
