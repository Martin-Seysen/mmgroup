import os
import sys
import tokenize
from pathlib import Path


def is_triple_quoted(string_token):
    """
    Return True if the token is a triple-quoted string.
    """
    # Remove any string prefixes like r, u, f, b, fr, etc.
    i = 0
    while i < len(string_token) and string_token[i].lower() in "rubf":
        i += 1

    quote_part = string_token[i:]
    return (
        quote_part.startswith('"""') or
        quote_part.startswith("'''")
    )


def is_raw_string(string_token):
    """
    Return True if the string has an r or R prefix.
    """
    prefix = ""
    i = 0
    while i < len(string_token) and string_token[i].lower() in "rubf":
        prefix += string_token[i].lower()
        i += 1
    return 'r' in prefix


def find_problem_strings(root_directory):
    root = Path(root_directory)

    for py_file in root.rglob("*.py"):
        try:
            with py_file.open("rb") as f:
                tokens = tokenize.tokenize(f.readline)
                for token in tokens:
                    if token.type == tokenize.STRING:
                        string_value = token.string

                        if (
                            is_triple_quoted(string_value)
                            and not is_raw_string(string_value)
                            and "\\" in string_value
                        ):
                            print(f"{py_file}:{token.start[0]}")
                            print("  Problematic string:")
                            print(f"  {string_value}\n")

        except Exception as e:
            print(f"Error reading {py_file}: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python find_strings.py <directory>")
        sys.exit(1)

    find_problem_strings(sys.argv[1])

