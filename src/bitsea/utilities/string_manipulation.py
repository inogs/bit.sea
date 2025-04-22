from __future__ import annotations

import re
from collections.abc import Iterable
from typing import List
from typing import Union


class ComposedString:
    """
    A `ComposedString` represents a string constructed from multiple substrings,
    including both plain strings and other `ComposedString` objects.

    This class simplifies the manipulation of strings containing parentheses.
    For example, the string:
        `((x + 1) * (y + 2)) + 3`
    can be represented as a `ComposedString` object made up of two substrings:
        `(x + 1) * (y + 2)` and ` + 3`.
    Furthermore, the first substring can itself be split further into three
    substrings:
        `x + 1`, ` * `, and `y + 2`.

    The `split_on_parenthesis` method can be used to convert a regular Python
    string into a hierarchical tree structure of `ComposedString` objects.
    Alternatively, you can directly create a `ComposedString` by providing a
    list of substrings, which are stored in the `elements` attribute.

    When a `ComposedString` is printed, it is converted into a regular string
    by concatenating its elements. Plain strings within the `ComposedString`
    are concatenated as is, while nested `ComposedString` objects are enclosed
    in parentheses.
    """

    def __init__(self, elements: Iterable[Union[str, ComposedString]]):
        self.elements = tuple(elements)

    def __str__(self) -> str:
        output = ""
        for element in self.elements:
            if isinstance(element, ComposedString):
                output += "(" + str(element) + ")"
            else:
                output += element
        return output

    def __repr__(self) -> str:
        return f"ComposedString({repr(self.elements)})"

    def __eq__(self, other):
        if not hasattr(other, "elements"):
            return NotImplemented
        return self.elements == other.elements

    @staticmethod
    def split_on_parenthesis(txt_str: str) -> ComposedString:
        """
        Given a common string, returns a new ComposedString that recursively
        splits it into a sequence of ComposedString objects around the
        parenthesis.
        """
        # A list to store the elements of the resulting ComposedString object
        pieces = []

        # A buffer to accumulate characters for the current substring being processed
        current_piece = ""

        # Tracks the number of currently open parentheses to determine nesting level
        current_level = 0

        open_parenthesis = "("
        close_parenthesis = ")"

        for char in txt_str:
            # Handle the occurrence of an open parenthesis
            if char == open_parenthesis:
                # If not already inside parentheses, save the current substring
                # as an independent piece and begin a new substring for the content inside
                # the parentheses.
                if current_level == 0:
                    if len(current_piece) > 0:
                        pieces.append(current_piece)
                        current_piece = ""
                # Otherwise, append the character to the current substring being processed
                else:
                    current_piece += char

                # Increase the nesting level to account for this new open parenthesis
                current_level += 1
                continue

            if char == close_parenthesis:
                # Decrease the nesting level to account for this closed parenthesis
                current_level -= 1

                # If the nesting level becomes negative, it means there are more
                # closing parentheses than opening ones, indicating an unbalanced string
                if current_level < 0:
                    raise ValueError(
                        f"Unbalanced parenthesis in string {txt_str}"
                    )

                # If exiting the outermost set of parentheses, convert the
                # accumulated substring into a new ComposedString and reset
                # the buffer
                if current_level == 0:
                    pieces.append(
                        ComposedString.split_on_parenthesis(current_piece)
                    )
                    current_piece = ""
                else:
                    current_piece += char
                continue

            current_piece += char

        # If the loop ends with unmatched open parentheses, the string is unbalanced
        if current_level != 0:
            raise ValueError("Unbalanced parenthesis in string: " + txt_str)

        # Append any remaining substring as the final piece if the string
        # did not end with a closing parenthesis
        if len(current_piece) > 0:
            pieces.append(current_piece)

        return ComposedString(pieces)


def _dequote(txt: str):
    """
    Extract the contents of a quoted field in a CSV file.

    Fields in CSV files may require quoting if they contain special characters
    like commas, quotes, or newlines. To escape quotes within the field, they
    are doubled. This function restores the original content by replacing each
    pair of quotes with a single quote.

    Note:
        This is a helper function used internally by `read_csv_txt` and should
        not be called directly.
    """
    # Check if an odd number of quotes appears in the field; if so, raise an
    # error
    odd_quotes = re.compile(r'(^|[^"])"("")*([^"]|$)')
    if odd_quotes.search(txt) is not None:
        raise ValueError(
            f"Invalid CSV file: The field {txt} contains an odd number of "
            f"quotes."
        )

    # This regular expression verifies that there is an even number of quotes
    # starting from the current position
    even_quotes = re.compile(r'("")+([^"]|$)')

    char_index = 0
    unquoted_field = ""
    while char_index < len(txt):
        # If the current character is not a quote, append it directly to the
        # result
        if txt[char_index] != '"':
            unquoted_field += txt[char_index]
            char_index += 1
            continue

        # Otherwise, count the number of consecutive quotes by matching the
        # pattern
        quotes = even_quotes.match(txt[char_index:])
        if quotes is None:
            raise ValueError(
                f"Invalid CSV file: The field {txt} contains an odd number of "
                f"quotes."
            )

        # The `group(0)` of the regular expression contains all consecutive
        # quotes, excluding the last non-quote character (if any). Count the
        # quotes explicitly.
        n_quotes = quotes.group(0).count('"')

        # Append half the number of quotes to the result and skip the
        # corresponding number of characters
        unquoted_field += '"' * (n_quotes // 2)
        char_index += n_quotes

    return unquoted_field


def read_csv_txt(txt: str) -> List[List[str]]:
    """
    Parse the content of a CSV file into a list of lists.

    Although Python's standard library includes a CSV reader, it cannot handle
    fields that exceed a specific length. While this is usually sufficient for
    most files, certain files (e.g., those containing WKT polygons) require the
    ability to process arbitrarily long fields. This function addresses that
    limitation with a pure Python implementation of a CSV reader.

    Note:
        This implementation is less efficient and less robust than the standard
        library. Use the default CSV reader when possible for better performance
        and reliability.
    """
    lines = [[]]
    current_line = lines[-1]

    # This regular expression identifies a quoted field. In such cases, quotes
    # within the field are doubled. The regular expression ensures that the
    # field contains either characters other than quotes or pairs of quotes.
    # It expects the quoted field to end with a single quote followed by a
    # comma, newline, or the end of the file.
    quoted_field = re.compile(r'^"(?P<field>([^"]|(""))*)"(?P<end>([,\n]|$))')

    # For unquoted fields, identify the end of the field by finding the first
    # comma or newline. If neither is found, the field ends at the end of the
    # file.
    unquoted_field = re.compile(r'^(?P<field>[^",\n]*?)(?P<end>([,\n]|$))')

    # The file_end pattern checks if we have reached the end of the file,
    # accommodating possible empty lines at the end.
    file_end = re.compile(r"(\s|\n)*$")

    while True:
        # Check the first character of the current record
        first_char = txt[0]

        # If it is a quote, we have found a quoted record
        is_a_quoted_field = first_char == '"'

        # Choose an appropriate regular expression
        field_re = quoted_field if is_a_quoted_field else unquoted_field

        # and apply it on the field
        field_match = field_re.match(txt)
        if field_match is None:
            raise ValueError(f'Invalid CSV file; error while reading "{txt}"')

        # This is the content of the field (beside the quotes)
        field_content = field_match.group("field")

        # and this is the amount of characters of the file we have read to
        # reach the end of the record (comma or newline at the end included).
        read_chars = len(field_match.group(0))

        # Remove from the file the characters we have read
        txt = txt[read_chars:]

        if is_a_quoted_field:
            field_content = _dequote(field_content)
        current_line.append(field_content)

        field_end = field_match.group("end")

        # If we reach the end of the file, we stop the loop
        if file_end.match(txt):
            # If we reached the end of the file but the last char was a comma,
            # there is a last empty field we have to save
            if field_end == ",":
                current_line.append("")
            break

        # If the field_end was a newline, we need to prepare a new line
        # where to save the content of the next field. This must be done after
        # the break for the end of file, otherwise we may prepare an empty line
        if field_end != ",":
            lines.append([])
            current_line = lines[-1]

    # Check that all the lines have the same length
    line_lengths = sorted(list(set(map(len, lines))))
    if len(line_lengths) > 1:
        raise ValueError(
            "Find lines of different length inside the CSV file: {}".format(
                ", ".join([str(line_length) for line_length in line_lengths])
            )
        )

    return lines
