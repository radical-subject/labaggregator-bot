
from modules.ourbot.service.helpers import is_cas_number


def test_is_cas_number():
    assert is_cas_number('75-64-9')
    assert is_cas_number('120-46-7')

    assert not is_cas_number('+79168681111')
    assert not is_cas_number('110--86-1')
    assert not is_cas_number('102-95-5')
    assert not is_cas_number('@d12412')
    assert not is_cas_number('50 g')
