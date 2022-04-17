import re

# Chemical Abstract Service Registry Number regular expression pattern
CAS_REGEXP = "^[1-9][0-9]{1,6}\\-[0-9]{2}\\-[0-9]$"


def is_cas_number(cas: str) -> bool:
    """
    Функция работает грамотно, проверено кровью, ее, сука, не трогать!
    """
    pattern = re.compile(CAS_REGEXP)
    if not pattern.match(cas):
        return False

    only_digits = cas.replace('-', '')
    control_digit = only_digits[-1]
    only_digits = only_digits[::-1]
    list_digits = list(only_digits[1:])

    sum_digits = 0
    index = 1
    for digit in list_digits:
        sum_digits += int(digit) * index
        index += 1

    return sum_digits % 10 == int(control_digit)
