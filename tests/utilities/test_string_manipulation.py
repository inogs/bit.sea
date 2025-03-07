from bitsea.utilities.string_manipulation import ComposedString
from bitsea.utilities.string_manipulation import read_csv_txt


def test_composed_string_init1():
    test_case = "str0 (str1)((str2)+(str3))"
    c1 = ComposedString.split_on_parenthesis(test_case)
    assert len(c1.elements) == 3
    assert test_case == str(c1)


def test_composed_string_init2():
    test_case = "str0 (str1) str3"
    c1 = ComposedString.split_on_parenthesis(test_case)
    assert len(c1.elements) == 3
    assert test_case == str(c1)


def test_composed_string_no_parenthesis():
    test_case = "str1 + str2"
    c1 = ComposedString.split_on_parenthesis(test_case)
    assert test_case == str(c1)


def test_composed_string_all_parenthesis():
    test_case = "str1 + str2"
    c1 = ComposedString.split_on_parenthesis(test_case)
    c2 = ComposedString.split_on_parenthesis("(" + test_case + ")")

    assert len(c1.elements) == 1
    assert len(c2.elements) == 1

    assert c2.elements[0] == c1


def test_read_csv_txt1():
    test_case = "str1,str2,str3\nstr4,str5,str6\n"
    c1 = read_csv_txt(test_case)
    assert len(c1) == 2
    assert c1[0] == ["str1", "str2", "str3"]
    assert c1[1] == ["str4", "str5", "str6"]


def test_read_csv_txt2():
    # Check that it works also if the last char is not a newline
    test_case = "str1,str2,str3\nstr4,str5,str6"
    c1 = read_csv_txt(test_case)
    assert len(c1) == 2
    assert c1[0] == ["str1", "str2", "str3"]
    assert c1[1] == ["str4", "str5", "str6"]


def test_read_csv_txt3():
    # Check that it works also if one field misses at the end of the line
    test_case = "str1,str2,\nstr4,str5,str6"
    c1 = read_csv_txt(test_case)
    assert len(c1) == 2
    assert c1[0] == ["str1", "str2", ""]
    assert c1[1] == ["str4", "str5", "str6"]


def test_read_csv_txt4():
    # Check that it works also if one field misses in the middle of the line
    test_case = "str1,,str3\nstr4,str5,str6"
    c1 = read_csv_txt(test_case)
    assert len(c1) == 2
    assert c1[0] == ["str1", "", "str3"]
    assert c1[1] == ["str4", "str5", "str6"]


def test_read_csv_txt5():
    # Check that it works also if one field misses at the end of the file
    test_case = "str1,str2,str3\nstr4,str5,\n"
    c1 = read_csv_txt(test_case)
    assert len(c1) == 2
    assert c1[0] == ["str1", "str2", "str3"]
    assert c1[1] == ["str4", "str5", ""]
