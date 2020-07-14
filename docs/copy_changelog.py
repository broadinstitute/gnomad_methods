#!/usr/bin/env python3

"""
Copy CHANGELOG.md from repository root and remove "Unreleased" section.
"""

import os


docs_directory = os.path.dirname(__file__)

changelog_path = os.path.join(docs_directory, "..", "CHANGELOG.md")
docs_changelog_path = os.path.join(docs_directory, "changelog.md")

with open(changelog_path) as f_in, open(docs_changelog_path, "w") as f_out:
    skip = False
    for line in f_in:
        is_heading = line.startswith("## ")
        if is_heading:
            skip = "unreleased" in line.lower()

        if not skip:
            f_out.write(line)
