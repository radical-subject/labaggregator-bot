
from modules.ourbot.service.helpers import is_CAS_number


def test_is_CAS_number():
    assert is_CAS_number('75-64-9')
    assert is_CAS_number('120-46-7')

    assert not is_CAS_number('+79168681111')
    assert not is_CAS_number('110--86-1')
    assert not is_CAS_number('102-95-5')
    assert not is_CAS_number('@d12412')
    assert not is_CAS_number('50 g')
