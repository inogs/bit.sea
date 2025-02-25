from bitsea.utilities.string_manipulation import ComposedString


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
