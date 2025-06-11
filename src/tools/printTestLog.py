#! /usr/bin/env python3

import re
from typing import Required
from rich.console import Console
from rich.repr import T
from rich.table import Table
from markdown import markdown
from bs4 import BeautifulSoup

def process_and_print_results(file_path):
    # Initialize console for rich output
    console = Console()

    # Read and parse the .md file
    with open(file_path, "r", encoding="utf-8") as file:
        md_content = file.read()

    # Extract the table section from the markdown
    table_pattern = re.compile(r"\|.*\|\n\|.*\|\n((?:\|.*\|\n)+)")
    table_match = table_pattern.search(md_content)
    if not table_match:
        raise ValueError("No table found in the markdown file.")

    # Extract header and rows
    lines = md_content.strip().split("\n")
    header_line_index = next(i for i, line in enumerate(lines) if line.startswith("|"))  # Locate header
    header = lines[header_line_index].strip("|").split("|")
    rows = lines[header_line_index + 2:]  # Skip the header and separator line

    # Clean and extract colored content
    def clean_html(cell):
        soup = BeautifulSoup(markdown(cell), "html.parser")
        return soup.get_text()

    data = [row.strip("|").split("|") for row in rows]
    data = [[clean_html(cell.strip()) for cell in row] for row in data]

    # Create a Rich table
    table = Table(title="Test Results", show_lines=True)

    # Define table columns
    for col in header:
        table.add_column(col.strip())

    print(data)

    # Add rows to the table with conditional coloring
    for row in data:
        result = row[2].strip()
        result_color = "green" if "Success" in result else "red"
        result_text = f"[{result_color}]{result}[/{result_color}]"

        comp_time_dif = row[5].strip()
        if "Success" in result:
            if "-" in comp_time_dif:
                dif_color = "green"
            else:
                dif_color = "red"
        else:
            dif_color = "red"
        comp_time_dif_text = f"[{dif_color}]{comp_time_dif}[/{dif_color}]"

        tests = row[6].strip()

        test_result = row[7].strip()
        if "Success" in test_result:
            test_result_color = "green"
            test_result_text = f"[{test_result_color}]{test_result}[/{test_result_color}]"
        elif "Failed" in test_result:
            test_result_color = "red"
            test_result_text = f"[{test_result_color}]{test_result}[/{test_result_color}]"
        else :
            test_result_text = f"{test_result}"

        # Add the formatted row to the table
        table.add_row(
            row[0].strip(),
            row[1].strip(),
            result_text,
            row[3].strip(),
            row[4].strip(),
            comp_time_dif_text,
            tests,
            test_result_text,
        )

    # Render the table
    console.print(table)

if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser(description="Arguments for script printing results of test")
    argparser.add_argument("--file", type=str, help="file with results", required=True)

    args = argparser.parse_args()
    process_and_print_results( args.file )
